#cython: language_level=3
import numpy as np
cimport numpy as np
import re
import copy
import json
from bookend.core.cython_utils._element_graph import ElementGraph
import bookend.core.cython_utils._rnaseq_utils as ru # RNAseqMapping, ELdata, range_of_reads, get_gaps, get_source_dict, build_depth_matrix
from collections import deque, Counter
import cython
import time

cdef class EndRange:
    cdef public int endtype
    cdef public int left, right, peak, terminal, strand
    cdef public float weight
    cdef public str tag
    """Represents a reference point for a Start or End site."""
    def __init__(self, left, right, peak, weight, endtype):
        self.left, self.right, self.peak, self.weight, self.endtype = left, right, peak, weight, endtype
        self.tag, self.strand = [('S', 1), ('E', 1), ('S', -1), ('E', -1)][self.endtype]
        if self.endtype in [0, 3]:
            self.terminal = left
        else:
            self.terminal = right
    
    def __repr__(self):
        strand = ['.','+','-'][self.strand]
        return '{}{}{} ({})'.format(self.tag, strand, self.peak, self.weight)
    
    def write_as_bed(self, chrom):
        print('{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, self.left, self.right, self.tag, self.weight, {-1:'-',1:'+',0:'.'}[self.strand]))

cdef class Locus:
    cdef public int chrom, leftmost, rightmost, extend, end_extend, number_of_elements, min_overhang, chunk_number, oligo_len
    cdef public bint naive,  use_attributes, ignore_ends
    cdef public tuple reads, frags
    cdef public float weight, bases, raw_bases, minimum_proportion, cap_bonus, intron_filter
    cdef public dict J_plus, J_minus, end_ranges, source_lookup
    cdef public set branchpoints
    cdef public list transcripts, traceback, sources
    cdef public object graph
    cdef EndRange nullRange
    cdef public np.ndarray depth_matrix, cov_plus, cov_minus, depth, read_lengths, member_lengths, frag_len, frag_by_pos, strand_array, weight_array, rep_array, membership, overlap, information_content, member_content, frag_strand_ratios, member_weights
    def __init__(self, chrom, chunk_number, list_of_reads, extend=50, end_extend=100, min_overhang=3, reduce=True, minimum_proportion=0.01, cap_bonus=5, complete=False, verbose=False, naive=True, intron_filter=0.15, use_attributes=False, oligo_len=20, ignore_ends=False):
        self.nullRange = EndRange(-1, -1, -1, -1, -1)
        self.oligo_len = oligo_len
        self.transcripts = []
        self.traceback = []
        self.branchpoints = set()
        self.member_weights = np.empty(0)
        self.chunk_number = chunk_number
        self.naive = naive
        self.minimum_proportion = minimum_proportion
        self.intron_filter = intron_filter
        self.min_overhang = min_overhang
        self.chrom = chrom
        self.cap_bonus = cap_bonus
        self.use_attributes = use_attributes
        self.ignore_ends = ignore_ends
        if len(list_of_reads) > 0:
            self.leftmost, self.rightmost = ru.range_of_reads(list_of_reads)
            if self.ignore_ends:
                for read in list_of_reads:
                    read.s_tag, read.e_tag = False, False
                    if len(read.splice)==0:read.strand = 0
            
            self.reads = tuple(list_of_reads) # Cannot be mutated            
            self.read_lengths = np.array([r.get_length()+self.oligo_len*r.s_tag+self.oligo_len*r.e_tag for r in self.reads], dtype=np.int32)
            self.raw_bases = np.sum(self.read_lengths * np.array([read.weight for read in self.reads]))
            self.sources = ru.get_sources(list_of_reads)
            self.source_lookup = ru.get_source_dict(self.sources)
            self.extend = extend
            self.end_extend = end_extend
            self.depth_matrix, self.J_plus, self.J_minus = ru.build_depth_matrix(self.leftmost, self.rightmost, self.reads, self.cap_bonus, self.use_attributes)
            self.prune_junctions()
            self.generate_branchpoints()
            if type(self) is AnnotationLocus:
                self.traceback = [set([i]) for i in range(len(self.reads))]
                empty = self.build_membership_matrix(0)
                if not empty:
                    self.filter_by_reps(self.minreps)
                    if self.membership.shape[0] > 0:
                        self.build_overlap_matrix(reduce=False, ignore_ends=True)
            else:
                empty = self.build_membership_matrix()
                if not empty:
                    self.build_overlap_matrix(reduce)
                    self.build_graph(reduce)
    
    def __len__(self):
        return self.rightmost - self.leftmost
    
    def __add__(self, other):
        return self.raw_bases + other.raw_bases
    
    def __repr__(self):
        symbols = {-1:'-', 0:' ', 1:'+', 2:'^'}
        summary_string = '<{} ({})>\n'.format(str(type(self)).split("'")[-2], self.number_of_elements)
        for l in range(self.number_of_elements):
            members = ''.join([symbols[i] for i in self.membership[l,:]])
            overlap = ''.join([symbols[i] for i in self.overlap[l,:]])
            indices = str(l) + (' ' if l < 100 else '') + (' ' if l < 10 else '')
            summary_string += '{} |{}|\t|{}|\n'.format(indices, members, overlap)
        
        return summary_string
    
    cpdef void prune_junctions(self):
        """If a splice junction represents < minimum_proportion of junction-spanning reads
        that either share its donor or its acceptor site, remove it."""
        cdef dict jdict, leftsides, rightsides
        cdef list keys, spans
        cdef str k
        cdef int i,l,r
        cdef set prune
        cdef np.ndarray counts
        cdef float total
        cdef bint passes
        for jdict in [self.J_plus, self.J_minus]:
            if jdict:
                prune = set()
                keys = sorted(list(jdict.keys()))
                spans = [self.string_to_span(k) for k in keys]
                leftsides = {}
                rightsides = {}
                for i in range(len(keys)):
                    l,r = spans[i]
                    leftsides[l] = leftsides.get(l, [])+[keys[i]]
                    rightsides[r] = rightsides.get(r, [])+[keys[i]]
                
                for l in leftsides.keys():
                    if len(leftsides[l]) > 1:
                        counts = np.array([jdict[k] for k in leftsides[l]], dtype=np.float32)
                        total = np.sum(counts)
                        prune.update([k for k,fails in zip(leftsides[l], counts/total < self.minimum_proportion) if fails])
                    
                for r in rightsides.keys():
                    if len(rightsides[r]) > 1:
                        counts = np.array([jdict[k] for k in rightsides[r]], dtype=np.float32)
                        total = np.sum(counts)
                        prune.update([k for k,fails in zip(rightsides[r], counts/total < self.minimum_proportion) if fails])
                
                for k in prune:
                    del jdict[k]
    
    cpdef void generate_branchpoints(self):
        """Estimate strand-specific coverage, then use it to generate an ordered array of
        positions where branching could occur in the overlap graph."""
        cdef:
            np.ndarray strandratio, covstranded, strandedpositions, pos, vals, value_order
            Py_ssize_t Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn, covrow, i
            float threshold_depth, cumulative_depth, cutoff
            int l, r, p
            EndRange rng
            set prohibited_positions
            list gaps
        
        Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn = range(11)
        covstranded = np.sum(self.depth_matrix[(covp,covm),:],axis=0)
        strandedpositions = np.where(covstranded > 0)[0]
        if strandedpositions.shape[0] > 0: # Some reads to inform the strand
            strandratio = np.array(np.interp(range(covstranded.shape[0]), strandedpositions, self.depth_matrix[covp,strandedpositions]/covstranded[strandedpositions]),dtype=np.float32)
            self.cov_plus = self.depth_matrix[covp,:] + self.depth_matrix[covn,:]*strandratio
            self.cov_minus = self.depth_matrix[covm,:] + self.depth_matrix[covn,:]*(1-strandratio)
            self.depth = self.cov_plus + self.cov_minus
        else:
            self.cov_plus = self.depth_matrix[covn,:]*.5
            self.cov_minus = self.depth_matrix[covn,:]*.5
            self.depth = self.depth_matrix[covn,:]
        
        self.branchpoints = set()
        self.end_ranges = dict()
        prohibited_positions = set()
        for j in self.J_plus.keys():
            self.branchpoints.update(list(self.string_to_span(j)))
        
        for j in self.J_minus.keys():
            self.branchpoints.update(list(self.string_to_span(j)))
        
        if self.min_overhang > 0:
            for p in self.branchpoints:
                prohibited_positions.update(range(p-self.min_overhang, p+self.min_overhang+1))
        
        for endtype in [Sp, Ep, Sm, Em]:
            pos = np.where(self.depth_matrix[endtype,]>0)[0]
            if endtype in [Sp, Ep]:
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_plus[pos]
            else:
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_minus[pos]
            
            self.end_ranges[endtype] = self.make_end_ranges(pos, vals, endtype)
            self.end_ranges[endtype] = [rng for rng in self.end_ranges[endtype] if rng.peak not in prohibited_positions]
            self.branchpoints.update([rng.terminal for rng in self.end_ranges[endtype]])
        
        for endtype in [Sp, Ep, Sm, Em]:
            for rng in self.end_ranges[endtype]:
                prohibited_positions.update(range(rng.left, rng.right+1))
        
        self.add_gaps(prohibited_positions)
        if 0 not in prohibited_positions:self.branchpoints.add(0)
        if len(self) not in prohibited_positions:self.branchpoints.add(len(self))
    
    cpdef void add_gaps(self, set prohibited_positions):
        """Updates branchpoints to include the starts and ends of coverage gaps"""
        cdef float cutoff
        cdef list gaps_plus, gaps_minus
        cdef (int, int) block, maxdelta
        if np.sum(self.cov_plus)>0:
            cutoff = max(1., np.mean(self.cov_plus[self.cov_plus>0])*self.minimum_proportion)
            gaps_plus = ru.get_gaps(self.cov_plus, self.extend, cutoff)
            for block in gaps_plus: # Iterate over plus-stranded gaps
                l, r = block
                if l not in prohibited_positions: # Potential unlabeled e_plus
                    self.branchpoints.add(l)
                
                if r not in prohibited_positions: # Potential unlabeled s_plus
                    self.branchpoints.add(r)
        
        if np.sum(self.cov_minus)>0:
            cutoff = max(1., np.mean(self.cov_minus[self.cov_minus>0])*self.minimum_proportion)
            gaps_minus = ru.get_gaps(self.cov_minus, self.extend, cutoff)
            for block in gaps_minus: # Iterate over minus-stranded gaps
                l, r = block
                if l not in prohibited_positions:
                    self.branchpoints.add(l)
                
                if r not in prohibited_positions:
                    self.branchpoints.add(r)
    
    cpdef list make_end_ranges(self, np.ndarray pos, np.ndarray vals, int endtype):
        """Returns a list of tuples that (1) filters low-signal positions
        and (2) clusters high-signal positions within self.end_extend.
        Returns list of (l,r) tuples demarking the edges of clustered end positions."""
        cdef:
            np.ndarray value_order
            EndRange e
            int p, maxp
            float cumulative, threshold, v, maxv, weight
            list filtered_pos, end_ranges
            (int, int) current_range
        
        if len(pos) == 0:
            return []
        elif len(pos) == 1:
            return [EndRange(pos[0], pos[0]+1, pos[0], vals[0], endtype)]
        
        passes_threshold = vals > self.minimum_proportion*np.sum(vals)
        # filtered_pos = sorted([(p,v) for p,v in zip(pos[value_order[:i]], vals[value_order[:i]])])
        filtered_pos = [(p,v) for p,v in zip(pos[passes_threshold], vals[passes_threshold])]
        p,v = filtered_pos[0]
        maxv = v
        maxp = p
        weight = v
        current_range = (p, p+1)
        end_ranges = []
        for p,v in filtered_pos[1:]:
            if p - self.end_extend <= current_range[1]:
                current_range = (current_range[0], p+1)
                weight += v
                if v > maxv or (v == maxv and endtype in [1,2]):
                    maxv, maxp = v, p
            else:
                e = EndRange(current_range[0], current_range[1], maxp, weight, endtype)
                end_ranges.append(e)
                current_range = (p, p+1)
                maxp = p
                maxv = v
                weight = v
        
        e = EndRange(current_range[0], current_range[1], maxp, weight, endtype)
        end_ranges.append(e)
        return end_ranges
    
    cdef int end_of_cluster(self, int pos, list end_ranges):
        """Returns the terminal position of the EndRange object that
        contains pos, if one exists. Else returns -1"""
        cdef EndRange rng
        for rng in end_ranges:
            if pos >= rng.left and pos <= rng.right:
                return rng.terminal
        
        return -1
    
    cdef EndRange get_end_cluster(self, int pos, list end_ranges):
        """Returns the most common position of the EndRange object that
        contains pos, if one exists. Else returns -1"""
        cdef EndRange rng
        for rng in end_ranges:
            if pos >= rng.left and pos <= rng.right:
                return rng
        
        return self.nullRange
    
    cdef str span_to_string(self, (int, int) span):
        """Converts a tuple of two ints to a string connected by ':'"""
        return '{}:{}'.format(span[0], span[1])
    
    cdef (int, int) string_to_span(self, str string):
        """Converts a string from span_to_string() back into a span"""
        cdef list splitstring = string.split(':')
        return (int(splitstring[0]), int(splitstring[1]))
    
    cpdef bint build_membership_matrix(self, float threshold=1):
        """After branchpoints are identified, populate a table that stores information about
        each read's membership within each frag:
            (1) read overlaps with frag
            (-1) read is incompatible with this frag
            (0) frag membership not determined
        
        The full table is produced by one pass through the reads. During this pass
        values are produced for each nucleotide of the locus object as:
            self.frag_by_pos - the unique fragment that each nucleotide belongs to
        After this, the table is compressed to produce:
            self.reduced_membership - a table of unique membership patterns
            self.membership_weights - (1 per row of reduced_membership) number of reads
        """
        cdef Py_ssize_t a, b, i, j, number_of_reads, number_of_frags, source_plus, source_minus, sink_plus, sink_minus
        cdef int last_rfrag, l, r, tl, tr, l_overhang, r_overhang, lfrag, rfrag, pos, locus_length
        cdef np.ndarray membership, strand_array, weight_array, discard, lengths, keep
        cdef char s
        cdef list temp_frags, bp_positions
        cdef set discard_frags
        cdef (int, int) block, span
        cdef str junction_hash
        
        bp_positions = sorted(list(self.branchpoints))
        temp_frags = []
        for i in range(len(bp_positions)-1):
            block = (bp_positions[i], bp_positions[i+1])
            temp_frags.append(block)
        
        self.frags = tuple(temp_frags) # Immutable after init, makes a vertex between all pairs of nonterminal branchpoint positions
        self.frag_len = np.array([b-a for a,b in self.frags]+[self.oligo_len]*4, dtype=np.int32)
        self.frag_strand_ratios = np.full(self.frag_len.shape[0], -1, dtype=np.float32)
        self.frag_strand_ratios[-4:-2] = 1.
        self.frag_strand_ratios[-2:] = 0.
        for i in range(len(self.frags)):
            span = self.frags[i]
            fragcov_plus = np.sum(self.cov_plus[span[0]:span[1]])
            fragcov_minus = np.sum(self.cov_minus[span[0]:span[1]])
            if fragcov_plus + fragcov_minus > 0:
                self.frag_strand_ratios[i] = fragcov_plus/(fragcov_plus+fragcov_minus)
        
        # self.frag_len = np.array([b-a for a,b in self.frags]+[0,0,0,0], dtype=np.int32)
        self.frag_by_pos = np.full(shape=len(self), fill_value=-1, dtype=np.int32)
        for i in range(len(self.frags)):
            if i == 0:
                a = i
            else:
                a = self.frags[i][0]
            
            if i == len(self.frags) - 1:
                b = len(self.frag_by_pos)
            else:
                b = self.frags[i][1]
                
            self.frag_by_pos[a:b] = i
        
        number_of_reads = len(self.reads)
        number_of_frags = len(self.frags)
        membership = np.zeros((number_of_reads, number_of_frags+4), dtype=np.int8) # Container for membership of each vertex in each read (-1 False, 1 True, 0 Unmeasured)
        strand_array = np.zeros(number_of_reads, dtype=np.int8) # Container for strandedness of each read (-1 minus, 1 plus, 0 nonstranded)
        weight_array = np.zeros((number_of_reads,len(self.source_lookup)), dtype=np.float32)
        self.member_lengths = np.zeros(number_of_reads, dtype=np.int32)
        self.rep_array = np.ones(number_of_reads, dtype=np.int32)
        source_plus, sink_plus, source_minus, sink_minus = range(number_of_frags, number_of_frags+4)
        locus_length = len(self.frag_by_pos)
        cdef char [:, :] MEMBERSHIP = membership
        for i in range(number_of_reads): # Read through the reads once, cataloging frags and branchpoints present/absent
            last_rfrag = 0
            read = self.reads[i]
            s = read.strand
            strand_array[i] = s
            if s == 1:
                MEMBERSHIP[i, source_minus] = -1 # Read cannot have minus-stranded features
                MEMBERSHIP[i, sink_minus] = -1 # Read cannot have minus-stranded features
            elif s == -1:
                MEMBERSHIP[i, source_plus] = -1 # Read cannot have plus-stranded features
                MEMBERSHIP[i, sink_plus] = -1 # Read cannot have plus-stranded features
            
            for j in range(len(read.ranges)): # Run through each block range of read
                block = read.ranges[j]
                l = block[0] - self.leftmost
                r = block[1] - self.leftmost
                # Get membership information about the block
                lfrag = self.frag_by_pos[l]
                rfrag = self.frag_by_pos[r-1]
                if self.min_overhang > 0:
                    if self.frag_len[lfrag] > self.min_overhang and l+self.min_overhang < locus_length:
                        lfrag = self.frag_by_pos[l+self.min_overhang-1]
                    
                    if self.frag_len[rfrag] > self.min_overhang and r-self.min_overhang >= 0:
                        rfrag = self.frag_by_pos[r-self.min_overhang]
                
                if j == 0: # Starting block
                    if s == 1 and read.s_tag: # Left position is a 5' end
                        # Sp, Ep, Sm, Em
                        tl = self.end_of_cluster(l, self.end_ranges[0])
                        if tl >= 0: # 5' end is in an EndRange
                            MEMBERSHIP[i, source_plus] = 1 # Add s+ to the membership table
                            l = tl
                            lfrag = self.frag_by_pos[l]
                            MEMBERSHIP[i, 0:lfrag] = -1 # Read cannot extend beyond source
                        else:
                            self.read_lengths[i] -= self.oligo_len
                    elif s == -1 and read.e_tag: # Left position is a 3' end
                        tl = self.end_of_cluster(l, self.end_ranges[3])
                        if tl >= 0:
                            MEMBERSHIP[i, sink_minus] = 1 # Add t- to the membership table
                            l = tl
                            lfrag = self.frag_by_pos[l]
                            MEMBERSHIP[i, 0:lfrag] = -1 # Read cannot extend beyond sink
                        else:
                            self.read_lengths[i] -= self.oligo_len
                
                if j == len(read.ranges)-1: # Ending block
                    if s == 1 and read.e_tag: # Right position is a 3' end
                        tr = self.end_of_cluster(r, self.end_ranges[1])
                        if tr >= 0:
                            MEMBERSHIP[i, sink_plus] = 1 # Add t+ to the membership table
                            r = tr
                            rfrag = self.frag_by_pos[r-1]
                            MEMBERSHIP[i, (rfrag+1):number_of_frags] = -1 # Read cannot extend beyond sink
                        else:
                            self.read_lengths[i] -= self.oligo_len
                    elif s == -1 and read.s_tag: # Right position is a 5' end
                        tr = self.end_of_cluster(r, self.end_ranges[2])
                        if tr >= 0:
                            MEMBERSHIP[i, source_minus] = 1 # Add s- to the membership table
                            r = tr
                            rfrag = self.frag_by_pos[r-1]
                            MEMBERSHIP[i, (rfrag+1):number_of_frags] = -1 # Read cannot extend beyond source
                        else:
                            self.read_lengths[i] -= self.oligo_len
                if lfrag > rfrag: # Reassignment of ends caused lfrag and rfrag to be out of order
                    if len(read.ranges) > 1:
                        if j == 0 and read.splice[j]: # The right border is a splice junction, structural violation
                            MEMBERSHIP[i,:] = -1 # Read is a violation, remove
                            break
                        elif j == len(read.ranges)-1 and read.splice[j-1]: # The left border is a splice junction, structural violation
                            MEMBERSHIP[i,:] = -1 # Read is a violation, remove
                            break

                    if (read.s_tag and s == 1) or (read.e_tag and s == -1): # The left border was updated
                        rfrag = lfrag
                    else: # The right border was updated
                        lfrag = rfrag
                
                MEMBERSHIP[i, lfrag:(rfrag+1)] = 1 # Add all covered frags to the membership table
                if j > 0: # Processing a downstream block
                    if read.splice[j-1]: # The gap between this block and the last was a splice junction
                        # Check that this junction is in the list of splice junctions
                        # If it was filtered out, this read is invalid and should be removed.
                        span = (self.frags[last_rfrag][1], self.frags[lfrag][0])
                        junction_hash = self.span_to_string(span)
                        if s == 1 and junction_hash not in self.J_plus:
                            MEMBERSHIP[i,:] = -1 # Read contains a filtered junction, remove
                            break
                        elif s == -1 and junction_hash not in self.J_minus:
                            MEMBERSHIP[i,:] = -1 # Read contains a filtered junction, remove
                            break
                        else:
                            MEMBERSHIP[i, (last_rfrag+1):lfrag] = -1 # All frags in the intron are incompatible
                
                last_rfrag = rfrag
            
            # Calculate the length (bases) of the element
            self.member_lengths[i] = np.sum(self.frag_len[membership[i,:] == 1])
            if self.member_lengths[i] > 0:
                weight_array[i, self.source_lookup[read.source]] += read.weight * self.read_lengths[i] / self.member_lengths[i]
        
        discard_frags = set()
        if threshold > 0:
            for i in range(len(self.frags)):
                l,r = self.frags[i]
                frag_depth = self.depth[l:r]
                if not passes_threshold(frag_depth, self.extend, threshold): # This frag has too large of a gap
                    discard_frags.add(i)
        
        discard = np.array(sorted(list(discard_frags)), dtype=np.int32)
        membership[np.sum(membership[:,discard]==1,axis=1) > 0,:] = -1 # Discard all elements with a discarded frag as a member
        if self.naive:
            weight_array = np.sum(weight_array,axis=1,keepdims=True)
        
        keep = np.sum(membership[:,:-4]==1,axis=1)>0
        self.membership = membership[keep,:]
        self.weight_array = weight_array[keep,:]
        self.strand_array = strand_array[keep]
        self.rep_array = self.rep_array[keep]
        self.member_lengths = self.member_lengths[keep]
        if not np.any(keep):
            return True
        
        self.reduce_membership()
        self.filter_members_by_strand()
        self.weight = np.sum(self.weight_array)
        self.number_of_elements = self.membership.shape[0]
        self.information_content = get_information_content(self.membership)
        self.member_content = get_member_content(self.membership)
        self.bases = np.sum(np.sum(self.weight_array, axis=1)*self.member_lengths)
        return False
    
    cdef void filter_members_by_strand(self):
        """If a read is in a region a region with >1-minimum_proportion coverage
        of a specific strand, assign this strand to the read. If the read is
        the opposite strand, discard it."""
        cdef float strand_ratio
        cdef np.ndarray lengths
        for i in range(self.membership.shape[0]):
            lengths = self.frag_len[self.membership[i,:] == 1]
            self.member_lengths[i] = np.sum(lengths)
            if self.member_lengths[i] > 0:
                strand_ratio = np.sum(self.frag_strand_ratios[self.membership[i,:]==1]*lengths)/self.member_lengths[i]
                if strand_ratio < self.minimum_proportion: # Minus-stranded region
                    if self.strand_array[i] == 0:
                        self.strand_array[i] = -1
                        self.membership[i, -4:-2] = -1 # Read cannot have plus-stranded features
                    elif self.strand_array[i] == 1: # Plus-stranded read in minus-stranded region; discard
                        self.membership[i,:] = -1
                elif strand_ratio > 1-self.minimum_proportion: # Plus-stranded region
                    if self.strand_array[i] == 0:
                        self.strand_array[i] = 1
                        self.membership[i, -2:] = -1 # Read cannot have minus-stranded features
                    elif self.strand_array[i] == -1:
                        self.membership[i,:] = -1
        
        self.reduce_membership()
    
    cpdef void denoise(self):
        """Checks the competitors of each read. If the read (plus compatible reads)
        accumulates to <minimum_proportion of sum(compatible + incompatible), 
        remove the read."""
        cdef np.ndarray keep, competitors, compatible, informative
        cdef float in_weight, out_weight, proportion
        cdef tuple compmembers
        keep = np.ones(self.membership.shape[0], dtype=np.bool)
        for index in range(self.membership.shape[0]):
            if keep[index]:
                competitors = self.get_competitors(index)
                if len(competitors) > 0:
                    compatible = self.get_compatible(index, competitors)
                    informative = np.where(np.apply_along_axis(np.all, 0, self.membership[compatible,:-4]==1))[0]
                    if len(informative) > 0:
                        # Subset for competitors that exclude the informative member(s) of the compatible set
                        competitors = competitors[np.sum(self.membership[competitors,:][:,informative]==-1, axis=1)>0]
                        # Subset for competitors that span the informative member(s)
                        compmembers = np.where(self.membership[competitors,:-4]==1)
                        competitors = competitors[sorted(set(compmembers[0][compmembers[1] < np.max(informative)]).intersection(set(compmembers[0][compmembers[1] > np.min(informative)])))]
                        in_weight = np.sum(self.weight_array[compatible,:])
                        out_weight = np.sum(self.weight_array[competitors,:])
                        proportion = in_weight / (in_weight+out_weight)
                        if proportion < self.minimum_proportion:
                            keep[compatible] = False
                        elif proportion > 1 - self.minimum_proportion:
                            keep[competitors] = False
        
        self.subset_elements(keep)
    
    cpdef np.ndarray get_competitors(self, int index):
        """Given an element index, return a list of all elements that 
        (1) share at least one member and 
        (2) are incompatible."""
        cdef np.ndarray incompatible, members, competitors
        incompatible =  np.where(self.overlap[:,index]==-1)[0]
        members = np.where(self.membership[index,:-4]==1)[0]
        competitors = incompatible[np.where(np.sum(self.membership[incompatible,:][:,members]==1, axis=1)>0)[0]]
        return competitors
        
    cpdef np.ndarray get_compatible(self, int index, np.ndarray competitors):
        cdef np.ndarray compatible
        compatible = np.where(np.logical_and(
            np.logical_or(
                self.overlap[:,index] > 0,
                self.overlap[index,:] > 0,
            ),
            np.all(self.overlap[:,competitors]<=0, axis=1)
        ))[0]
        return compatible
    
    cpdef void reduce_membership(self):
        """Given a matrix of membership values, 
        returns a [reduced_membership_matrix, weights] array
        such that all rows are unique and in sort order."""
        cdef np.ndarray reduced_membership, reverse_lookup, new_weights, new_strands, members_bool, new_lengths
        cdef list left_member, right_member, index, sort_triples, sorted_indices
        cdef (int, int, int) triple
        cdef Py_ssize_t i,v
        cdef bint member_weights_exists
        if self.membership.shape[0] > 1:
            reduced_membership, reverse_lookup = np.unique(self.membership, axis=0, return_inverse=True)
            new_weights = np.zeros(shape=(reduced_membership.shape[0], self.weight_array.shape[1]), dtype=np.float32)
            new_strands = np.zeros(shape=reduced_membership.shape[0], dtype=np.int8)
            new_reps = np.zeros(shape=reduced_membership.shape[0], dtype=np.int32)
            new_lengths = np.zeros(shape=reduced_membership.shape[0], dtype=np.int32)
            if not np.any(self.member_weights):
                member_weights_exists = False
            else:
                member_weights_exists = True
                new_member_weights = np.zeros(shape=(reduced_membership.shape[0],reduced_membership.shape[1]), dtype=np.float32)
            
            if type(self) is AnnotationLocus:
                new_traceback = []
                for i in range(reduced_membership.shape[0]):
                    new_traceback += [set()]
                
                for i,v in enumerate(reverse_lookup):
                    new_traceback[v].add(i)

            for i,v in enumerate(reverse_lookup):
                new_weights[v] += self.weight_array[i]
                new_strands[v] = self.strand_array[i]
                new_lengths[v] = self.member_lengths[i]
                new_reps[v] += self.rep_array[i]
                if member_weights_exists:
                    new_member_weights[v,:] += self.member_weights[i,:]
            
            members_bool = reduced_membership[:,[-4,-1]+list(range(0,reduced_membership.shape[1]-4))+[-3,-2]]==1
            number_of_members = np.sum(members_bool[:,2:-2],axis=1)
            left_member = np.argmax(members_bool, axis=1).tolist()
            right_member = (members_bool.shape[1]-1-np.argmax(members_bool[:,::-1], axis=1)).tolist()
            index = list(range(members_bool.shape[0]))
            sort_triples = sorted(list(zip(left_member, right_member, index)))
            sorted_indices = [triple[2] for triple in sort_triples if number_of_members[triple[2]] > 0]
            self.membership = reduced_membership[sorted_indices,:]
            self.weight_array = new_weights[sorted_indices,:]
            if not member_weights_exists:
                self.member_weights = np.full((self.membership.shape[0],self.membership.shape[1]), np.sum(self.weight_array,axis=1,keepdims=True))
                self.member_weights[self.membership==0] = 0
            else:
                self.member_weights = new_member_weights[sorted_indices,:]
            
            self.member_lengths = new_lengths[sorted_indices]
            self.strand_array = new_strands[sorted_indices]
            self.rep_array = new_reps[sorted_indices]
            if len(self.traceback) > 0:
                self.traceback = [new_traceback[i] for i in sorted_indices]
        
        if not np.any(self.member_weights): # member_weights still uninitialized
            self.member_weights = np.full((self.membership.shape[0],self.membership.shape[1]), np.sum(self.weight_array,axis=1,keepdims=True))
            self.member_weights[self.membership==0] = 0
    
    cpdef void filter_by_reps(self, int minreps=1):
        """Enforce that elements in the membership"""
        cdef np.ndarray keep
        if minreps > 1:
            keep = np.where(self.rep_array >= minreps)[0]
            self.subset_elements(keep)
    
    cpdef void build_overlap_matrix(self, bint reduce=True, bint ignore_ends=False):
        """Builds reduced-representation matrix where every read is (-1) incompatible,
        (1) overlapping, or (0) non-overlapping with every other read. Any read completely
        contained in one or more reads with more information are removed, with their
        weight distributed proportionally to the containers.
        """
        if ignore_ends:
            endless_matrix = remove_ends(self.membership)
            endless_info = get_information_content(endless_matrix)
            self.overlap = calculate_overlap(endless_matrix, endless_info, self.strand_array)
        else:
            self.overlap = calculate_overlap(self.membership, self.information_content, self.strand_array)
        
        if reduce:
            self.denoise()
            maxIC = self.membership.shape[1]
            self.resolve_containment()
            keep = np.where(np.sum(self.weight_array, axis=1) > 0)[0]
            self.subset_elements(keep)
    
    cpdef void subset_elements(self, np.ndarray keep):
        self.overlap = self.overlap[keep,:][:,keep]
        self.membership = self.membership[keep,:]
        self.weight_array = self.weight_array[keep,:]
        self.member_weights = self.member_weights[keep,:]
        self.rep_array = self.rep_array[keep]
        self.strand_array = self.strand_array[keep]
        self.information_content = self.information_content[keep]
        self.member_content = self.member_content[keep]
        self.member_lengths = self.member_lengths[keep]
        self.bases = np.sum(np.sum(self.weight_array, axis=1)*self.member_lengths)
        self.number_of_elements = self.membership.shape[0]
        if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
    
    cpdef void build_graph(self, reduce=True):
        """Constructs one or more graphs from 
        connection values (ones) in the overlap matrix.
        Additionally, stores the set of excluded edges for each node as an 'antigraph'
        """
        cdef np.ndarray updated
        if reduce: # Collapse linear chains prior to graph construction
            self.collapse_linear_chains()
        
        self.graph = ElementGraph(self.overlap, self.membership, self.weight_array, self.member_weights, self.strand_array, self.frag_len, self.naive)
    
    cpdef void assemble_transcripts(self, bint complete=False, bint collapse=True):
        cdef list reassigned_coverage
        cdef float total_coverage
        self.graph.assemble(self.minimum_proportion)
        if complete: # Discard all paths that aren't complete
            paths_to_remove = [i for i in range(len(self.graph.paths)) if not self.graph.paths[i].complete]
        else: # Still remove paths with internal gaps
            paths_to_remove = [i for i in range(len(self.graph.paths)) if self.graph.paths[i].has_gaps]
        
        self.graph.remove_paths(sorted(set(paths_to_remove)))
        self.graph.assign_weights()
        paths_to_remove = self.filter_truncations()
        if self.intron_filter > 0:
            paths_to_remove += self.filter_retained_introns()
        
        assigned_bases = np.array([path.bases for path in self.graph.paths])
        assigned_sum = np.sum(assigned_bases)
        updated_bases = np.array([path.bases for path in self.graph.paths])
        while np.sum(np.abs(updated_bases-assigned_bases)) < self.minimum_proportion * assigned_sum:
            assigned_bases = updated_bases
            self.graph.assign_weights()
            updated_bases = np.array([path.bases for path in self.graph.paths])
        
        for i in range(len(self.graph.paths)):
            path = self.graph.paths[i]
            self.transcripts.append(self.convert_path(path, i+1))
        
        self.add_transcript_attributes()
    
    cdef list filter_retained_introns(self):
        """Removes paths from the list of assembled paths if they
        contain at least one 'exonic' region that accumulates to
        < intron_filter (float 0-1) of the spliced + unspliced coverage."""
        cdef:
            list paths_to_remove
            char strand
            dict junctions_as_frags, junctions
            int i, number_of_paths, l, r, frag
            (int, int) junction, frag_pair
            str junction_hash, frag_hash
            float junction_count, intronic_coverage, spliced_coverage
            set members
        
        number_of_paths = len(self.graph.paths)
        paths_to_remove = []
        for strand in [-1, 1]:
            junctions_as_frags = {}
            junctions = self.J_plus if strand == 1 else self.J_minus
            for junction_hash, junction_count in junctions.items():
                # Convert all intron pair positions to frag indices
                block = self.string_to_span(junction_hash)
                l = block[0]
                r = block[1]
                frag_pair = (self.frag_by_pos[l], self.frag_by_pos[r-1])
                junctions_as_frags[self.span_to_string(frag_pair)] = junction_count      
            
            for i in range(number_of_paths):
                path = self.graph.paths[i]
                members = path.members
                if path.strand == strand:
                    for frag_hash in junctions_as_frags.keys():
                        # Check each like-stranded intron: is it fully intact?
                        block = self.string_to_span(frag_hash)
                        is_spliced_out = False
                        for frag in range(block[0],block[1]+1):
                            if frag not in members:
                                is_spliced_out = True
                        
                        if not is_spliced_out: # Intron is retained, apply filter
                            spliced_coverage = junctions_as_frags[frag_hash]
                            if path.coverage < self.intron_filter * (path.coverage + spliced_coverage):
                                paths_to_remove.append(i)
                                break
        
        return paths_to_remove
    
    cpdef list filter_truncations(self):
        """Removes paths from the list of assembled paths if they
        are a truncated version of a longer path in the set and they
        do not account for a majority of the reads of the set of transcripts
        compatible in this way."""
        cdef:
            list paths_to_remove, number_of_members, end_indices
            int index, n, f, i, j
            set members, nonmembers
            float container_reads
            np.ndarray paths_in_member_order
        
        paths_to_remove = []
        if len(self.graph.paths) <= 1:
            return paths_to_remove
        
        end_indices = sorted(list(self.graph.paths[0].end_indices))
        number_of_paths = len(self.graph.paths)
        # Visit each path once by increasing number of members
        number_of_members = [len(path.members) for path in self.graph.paths]
        paths_in_member_order = np.argsort(number_of_members)
        for i in range(number_of_paths-1):
            # Remove the start incompatibilities and check if the path
            # becomes contained in any other path.
            index = paths_in_member_order[i]
            path = self.graph.paths[index]
            members = copy.copy(path.members)
            nonmembers = copy.copy(path.nonmembers)
            if path.strand == 1:
                members.discard(end_indices[0]) # Remove S+
                members.discard(end_indices[1]) # Remove E+
            elif path.strand == -1:
                members.discard(end_indices[2]) # Remove S-
                members.discard(end_indices[3]) # Remove E-
            
            for f in range(path.LM): # Remove all incompatibilities upstream
                nonmembers.discard(f)
            
            for f in range(path.RM, end_indices[0]): # Remove all incompatibilities upstream
                nonmembers.discard(f)
            
            container_reads = 0
            for j in range(i+1, number_of_paths):
                other_index = paths_in_member_order[j]
                other_path = self.graph.paths[other_index]
                if other_path.strand == path.strand:
                    if len(members.difference(other_path.members)) == 0:
                        # All of path's members are in other_path
                        if len(nonmembers.intersection(other_path.members)) == 0:
                            if len(other_path.nonmembers.intersection(members)) == 0:
                                # path is fully contained in other_path
                                container_reads += np.sum(other_path.weights)
            
            # After all containers are found, compare path.reads to sum of all containers
            if np.sum(path.weights) < container_reads:
                paths_to_remove.append(index)

        return paths_to_remove
    
    cpdef void resolve_containment(self):
        """Given a overlap matrix, 'bubble up' the weight of
        all reads that have one or more 'contained by' relationships
        to other reads. Pass from highest complexity reads down, assigning
        weight proportional to the existing weight.
        The resulting matrix should contain only overlaps, exclusions, and unknowns."""
        cdef:
            np.ndarray containment, contained, IC_order, new_weights, containers, incompatible, weight_transform, compatibilities, informative, container_weights, container_proportions, total_container_weights, n_weights
            Py_ssize_t i
        
        containment = self.overlap==2 # Make a boolean matrix of which reads are contained in other reads
        np.put(containment, range(0,containment.shape[0]**2,containment.shape[0]+1), False, mode='wrap') # Blank out the diagonal (self-containments)
        contained = np.where(np.sum(containment, axis=1) > 0)[0] # Identify reads that are contained
        IC_order = contained[np.lexsort((-self.information_content[contained], -self.member_content[contained]))] # Rank them by decreasing number of members
        new_weights = np.copy(self.weight_array)
        for i in IC_order:
            containers = np.where(containment[i,:])[0]
            containers = containers[np.sum(new_weights[containers,:],axis=1)>0]
            # Get the set of reads incompatible with all containers but that do not exclude i
            incompatible =  np.where(np.logical_and(
                np.all(self.overlap[:,containers]==-1, axis=1),
                self.overlap[:,i] > 0
            ))[0]
            incompatible = incompatible[np.sum(new_weights[incompatible,:],axis=1)>0]
            if len(incompatible) == 0: # Special case, all weight goes to containers
                if len(containers) == 1:
                    new_weights[containers,:] += new_weights[i,:] * self.member_lengths[i]/self.member_lengths[containers]
                else: # Evaluate how much weight goes to each container
                    compatible = np.where(np.logical_and(
                        np.logical_or(
                            np.any(self.overlap[:,containers] > 0, axis=1),
                            np.any(self.overlap[containers,:] > 0, axis=0)
                        ),
                        np.all(self.overlap[:,incompatible]==-1, axis=1)
                    ))[0]
                    nonzero = np.where(new_weights[i,:] > 0)[0]
                    weight_to_add = np.zeros(shape=(len(containers),new_weights.shape[1]), dtype=np.float32)
                    weight_transform = self.member_lengths[i]/self.member_lengths[containers]
                    
                    compatibilities = self.overlap[:,containers][compatible,:] > -1
                    informative = np.where(np.sum(compatibilities,axis=1) < len(containers))[0]
                    container_weights = new_weights[containers,:]
                    for f in informative:
                        container_weights[compatibilities[f,:],:] += new_weights[compatible[f],:]
                        if compatible[f] in containers:
                            container_weights[containers==compatible[f],:] -= new_weights[containers[containers==compatible[f]],:]
                    
                    total_container_weights = np.sum(container_weights, axis=1)
                    container_proportions = total_container_weights / np.sum(total_container_weights)
                    for n in nonzero: # Each source that has reads of i is evaluated separately
                        n_weights = container_weights[:,n]
                        total = np.sum(n_weights)
                        weight = new_weights[i,n]
                        if total > 0:
                            weight_to_add[:,n] += weight * n_weights / total * weight_transform
                        else:
                            weight_to_add[:,n] += weight * container_proportions * weight_transform
                    
                    new_weights[containers,:] += weight_to_add
                
                new_weights[i,] = 0
                containment[:,i] = False
                containment[i,:] = False
    
    cpdef void add_transcript_attributes(self):
        """Populate the new read objects with diagnostic information
        to store in the GTF attributes column."""
        cdef:
            int first, last, s_pos, e_pos
            dict S_info, E_info
            list S_ranges, E_ranges
            EndRange S, E
        
        for T in self.transcripts:
            T.attributes['length'] = T.get_length()
            T.attributes['reads'] = round(T.weight, 2)
            S_info = {'S.reads':0, 'S.capped':0, 'S.left':0, 'S.right':0}
            E_info = {'E.reads':0, 'E.left':0, 'E.right':0}
            
            if T.strand != 0:
                if T.strand == 1:
                    s_pos = first = T.span[0]
                    e_pos = last = T.span[1]
                    S_ranges = self.end_ranges[0]
                    E_ranges = self.end_ranges[1]
                else:
                    s_pos = first = T.span[1]
                    e_pos = last = T.span[0]
                    S_ranges = self.end_ranges[2]
                    E_ranges = self.end_ranges[3]
                
                if T.s_tag:
                    S = self.get_end_cluster(first, S_ranges)
                    if S is not self.nullRange:
                        s_pos = S.peak
                        S_info['S.reads'] = round(S.weight, 2)
                        S_info['S.left'] = S.left
                        S_info['S.right'] = S.right
                        if s_pos != first: # S pos was replaced
                            if T.strand == 1:
                                T.ranges[0] = (s_pos, T.ranges[0][1])
                            else:
                                T.ranges[-1] = (T.ranges[-1][0], s_pos)
                
                if T.e_tag:
                    E = self.get_end_cluster(last, E_ranges)
                    if E is not self.nullRange:
                        e_pos = E.peak
                        E_info['E.reads'] = round(E.weight,2)
                        E_info['E.left'] = E.left
                        E_info['E.right'] = E.right
                        if e_pos != last: # S pos was replaced
                            if T.strand == 1:
                                T.ranges[-1] = (T.ranges[-1][0], e_pos)
                            else:
                                T.ranges[0] = (e_pos, T.ranges[0][1])
            
            T.attributes.update(S_info)
            T.attributes.update(E_info)
    
    cpdef void collapse_linear_chains(self):
        """Collapses chains of vertices connected with a single edge.
        """
        cdef int i, chain, parent
        cdef list keep = []
        cdef np.ndarray linear_chains, resolve_order
        cdef dict chain_parent = {}
        resolve_order = np.lexsort((-self.information_content, -self.member_content))
        linear_chains = find_linear_chains(self.overlap, resolve_order)
        for i,chain in enumerate(linear_chains):
            if chain == 0:
                keep.append(i)
            else:
                if chain in chain_parent.keys(): # chain exists, merge i into the parent
                    parent = chain_parent[chain]
                    self.merge_reads(i, parent)
                else: # chain doesn't yet exist, this is the parent
                    chain_parent[chain] = i
                    keep.append(i)
    
        if len(keep) < self.number_of_elements:
            self.subset_elements(np.array(sorted(keep)))
    
    cpdef void merge_reads(self, int child_index, int parent_index):
        """Combines the information of two read elements in the locus."""
        cdef char p, c, s
        cdef Py_ssize_t i
        cdef int combined_length
        for i in range(self.membership.shape[1]): # Iterate over columns of the membership table
            p = self.membership[parent_index,i]
            c = self.membership[child_index,i]  
            if p == 0: # If no info in parent, overwrite with child
                self.membership[parent_index, i] = c
            elif c != 0 and p != c: # Conflicting membership information
                raise Exception('Incompatible read pair: {}, {}'.format(child_index, parent_index))
        
        for i in range(self.overlap.shape[0]): # Iterate over columns/rows of the overlap matrix
            # rows
            p = self.overlap[parent_index, i]
            c = self.overlap[child_index, i]
            if p == 0 or c == -1:
                self.overlap[parent_index, i] = c
            
            # columns
            p = self.overlap[i, parent_index]
            c = self.overlap[i, child_index]
            if p == 0 or c == -1:
                self.overlap[i, parent_index] = c
        
        s = self.strand_array[parent_index]
        if s == 0:
            self.strand_array[parent_index] = self.strand_array[child_index]
        
        combined_length = np.sum(self.frag_len[self.membership[parent_index,:]==1])
        self.weight_array[parent_index,:] = (self.weight_array[child_index,:]*self.member_lengths[child_index] + self.weight_array[parent_index,:]*self.member_lengths[parent_index])/combined_length
        self.member_weights[parent_index,:] += self.member_weights[child_index,:]
        self.member_lengths[parent_index] = combined_length
        self.weight_array[child_index,:] = 0
    
    cpdef convert_path(self, element, transcript_number):
        """Prints a representation of an ElementGraph Element object
        by converting it first to an RNAseqMapping object."""
        cdef int chrom, source, N, m, n, l, r, last_member
        cdef char strand
        cdef bint s_tag, e_tag, capped, gap_is_splice
        cdef list ranges, splice, members
        cdef (int, int) frag, exon
        cdef set nonmembers
        cdef str gene_id, transcript_id, junction_hash
        cdef dict junctions
        cdef EndRange S, E
        gene_id = 'bookend.{}'.format(self.chunk_number)
        transcript_id = 'bookend.{}.{}'.format(self.chunk_number, transcript_number)
        members = sorted(element.members)
        nonmembers = element.nonmembers
        N = element.maxIC
        chrom = self.chrom
        source = 0
        strand = element.strand
        ranges = []
        splice = []
        l = r = last_member = -1
        gap_is_splice = True
        junctions = self.J_plus if strand == 1 else self.J_minus
        for m in members: # Convert membership into ranges
            # 0-indexed open doubles
            if m >= len(self.frags):
                break
            
            frag = self.frags[m]
            if last_member == -1: # Uninitialized
                l, r = frag
                # Update leftmost position if it matches an S/E branchpoint
                if element.s_tag and element.strand == 1:
                    S = self.get_end_cluster(l, self.end_ranges[0])
                    if S is not self.nullRange:
                        l = S.peak
                elif element.e_tag and element.strand == -1:
                    E = self.get_end_cluster(l, self.end_ranges[3])
                    if E is not self.nullRange:
                        l = E.peak
            elif last_member == m-1:
                r = frag[1]
            else: # A gap was jumped
                junction_hash = self.span_to_string((r, frag[0]))
                if junction_hash in junctions:
                    gap_is_splice = True
                else:
                    gap_is_splice = False
                
                splice.append(gap_is_splice)
                exon = (l+self.leftmost, r+self.leftmost)
                ranges.append(exon)
                l, r = frag
            
            last_member = m
        
        # Update rightmost position
        if element.s_tag and element.strand == -1:
            S = self.get_end_cluster(r, self.end_ranges[2])
            if S is not self.nullRange:
                r = S.peak+1
        elif element.e_tag and element.strand == 1:
            E = self.get_end_cluster(r, self.end_ranges[1])
            if E is not self.nullRange:
                r = E.peak+1
        
        exon = (l+self.leftmost, r+self.leftmost)
        ranges.append(exon)
        s_tag = element.s_tag
        e_tag = element.e_tag
        capped = False
        weight = np.sum(element.source_weights)
        elementAttributes = {}
        elementData = ru.ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)
        readObject = ru.RNAseqMapping(elementData, elementAttributes)
        readObject.attributes['gene_id'] = gene_id
        readObject.attributes['transcript_id'] = transcript_id
        readObject.attributes['cov'] = round(readObject.weight, 2)
        return readObject

##########################################

cdef class AnnotationLocus(Locus):
    """Processes sets of RNAseqMapping objects with behavior specific
    to 'bookend merge'"""
    cdef public bint ignore_reference_ends
    cdef public list ref_reads, high_confidence_transcripts
    cdef public int minreps, confidence
    cdef public float total_s, total_e, cap_percent
    cdef public AnnotationLocus ref_locus
    def __init__(self, chrom, chunk_number, list_of_reads, end_extend, min_overhang=0, minimum_proportion=0, cap_percent=0.1, intron_filter=0, ignore_reference_ends=True, minreps=1, confidence=-1):
        Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn = range(11)
        self.minreps = minreps
        self.confidence = confidence
        self.ignore_reference_ends = ignore_reference_ends
        self.cap_percent = cap_percent
        self.ref_reads, nonref_reads = self.get_annotation_info(list_of_reads)
        if len(nonref_reads) > 0:
            Locus.__init__(self, chrom, chunk_number, nonref_reads, 0, end_extend, min_overhang, True, minimum_proportion, 0, 1, False, False, True, intron_filter, False, False, True)
            self.cap_percent = cap_percent
            self.total_s = sum([er.weight for er in [self.end_ranges[Sp]]+[self.end_ranges[Sm]]])
            self.total_e = sum([er.weight for er in [self.end_ranges[Ep]]+[self.end_ranges[Em]]])
            self.transcripts = [self.transcript_from_membership(i) for i in range(self.membership.shape[0])]
            if self.confidence == -1:
                self.high_confidence_transcripts = []
            else:
                self.high_confidence_transcripts = list(np.where(self.rep_array >= self.confidence)[0])

    cpdef tuple get_annotation_info(self, list list_of_reads):
        """Updates features of the RNAseqMapping object with specific
        fields from its attribute dict, if they exist."""
        cdef list ref_reads, nonref_reads
        cdef float percent_capped, s_reads, c_reads
        ref_reads = []
        nonref_reads = []
        for read in list_of_reads:
            read.weight = float(read.attributes['TPM']) if 'TPM' in read.attributes.keys() else read.weight
            if read.is_reference:
                if self.ignore_reference_ends:
                    read.s_tag = read.e_tag = read.capped = False
                
                ref_reads.append(read) 
            else:
                s_reads = float(read.attributes.get('S.reads', 1))
                c_reads = float(read.attributes.get('S.capped', 1))
                percent_capped = c_reads/(s_reads+c_reads)
                read.capped = True if percent_capped >= self.cap_percent else False
                read.attributes['transcript_id'] = '{}.{}'.format(read.attributes['source'],read.attributes['transcript_id'])
                nonref_reads.append(read)
        
        return ref_reads, nonref_reads
    
    cpdef list identify_truncations(self):
        """Returns a list of member indices which are fully contained in
        longer transcript(s)."""
        containment = self.overlap==2 # Make a boolean matrix of which reads are contained in other reads
        return list(np.where(np.sum(containment, axis=1) > 1)[0])
    
    cpdef list identify_containers(self):
        """Returns a list of member indices which fully contain
        shorter transcript(s)."""
        containment = self.overlap==2 # Make a boolean matrix of which reads are contained in other reads
        return list(np.where(np.sum(containment, axis=0) > 1)[0])

    cpdef void filter_fused_and_truncated_annotations(self):
        """Identifies situations in the merged Membership matrix that
        represent putative false positive truncations (5' or 3' fragments)
        and false positive fusions (containment of two adjacent nonoverlapping genes).
        Requires that each must pass filtering criteria to be kept."""
        cdef:
            list truncation_indices, container_indices, to_remove, breaks
            int strand, i, c, pos
            set containers
            np.ndarray members, contained_segments, 
        
        if len(self.transcripts) == 0:
            return
        
        # TEST 1: Check if a truncation is capped (Short in Long with different TSS)
        # FILTER 1: Short must be capped
        truncation_indices = self.identify_truncations()
        to_remove = []
        if len(truncation_indices) == 0: # Nothing is contained in anything else
            return
        
        for i in truncation_indices:
            read = self.transcripts[i]
            strand = read.strand
            startpos = sorted([pos*strand for pos in read.span])[0]
            # At least 1 container has an upstream 5' end
            containers = set(np.where(self.overlap[i,:]==2)[0]).difference([i])
            upstream_start = any([sorted([pos*strand for pos in self.transcripts[c].span])[0] < startpos for c in containers])
            if upstream_start and not read.capped:
                to_remove.append(i)
        
        self.remove_transcripts(to_remove)
        # TEST 2: Check if a putative fusion exists (non-overlapping Short pair in Long)
        # FILTER 2: Long must have more weight than each Short
        to_remove = []
        container_indices = self.identify_containers()
        if len(container_indices) == 0: # Nothing is contained in anything else
            return
        
        container_indices = [container_indices[i] for i in np.argsort(-self.member_content[container_indices])]
        for i in container_indices: # Iterate over containers in decreasing length
            contained = list(set(np.where(self.overlap[:,i]==2)[0]).difference([i]))
            if len(contained) > 1:
                members = np.where(self.membership[i,:-4]==1)[0]
                contained_segments = self.membership[contained, :][:, members]
                breaks = find_breaks(contained_segments,False)
                if len(breaks) > 0: # i spans a split population of contained transcripts
                    c_weight = max([self.transcripts[c].weight for c in contained])
                    if self.transcripts[i].weight < c_weight:
                        to_remove.append(i)
        
        self.remove_transcripts(to_remove)
        # TEST 3: Check for putative fragments (spurious caps and end labels)
        # FILTER 3: Short must have more weight than each Long
        truncations = self.identify_truncations()
        to_remove = []
        if len(truncations) == 0: # Nothing is contained in anything else
            return
        
        for i in truncations:
            read = self.transcripts[i]
            containers = set(np.where(self.overlap[i,:]==2)[0]).difference([i])
            c_weight = max([self.transcripts[c].weight for c in containers])
            if read.weight < c_weight:
                to_remove.append(i)
        
        self.remove_transcripts(to_remove)
    
    cpdef remove_transcripts(self, list to_remove):
        """Reduce the locus to only the set of elements
        not present in the list to_remove."""
        cdef list keep
        cdef int r, k
        to_remove = [r for r in to_remove if r not in self.high_confidence_transcripts]
        if len(to_remove) == 0:
            return
        
        keep = [i for i in range(self.number_of_elements) if i not in to_remove]
        self.number_of_elements = len(keep)
        self.overlap = self.overlap[keep,:][:,keep]
        self.membership = self.membership[keep,:]
        self.weight_array = self.weight_array[keep,:]
        self.rep_array = self.rep_array[keep]
        self.strand_array = self.strand_array[keep]
        self.information_content = self.information_content[keep]
        self.member_content = self.member_content[keep]
        self.transcripts = [self.transcripts[k] for k in keep]
        self.bases = np.sum(np.sum(self.weight_array, axis=1)*self.member_lengths)
        if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
        if self.confidence >= 0: self.high_confidence_transcripts = list(np.where(self.rep_array >= self.confidence)[0])

    cpdef end_ratio(self, int index):
        """Calculates the average deviation from an equal ratio
        of Start and End read proportions across all reads of the meta-assembly."""
        traceback = self.traceback[index]
        average_ratio = 0
        for i in traceback:
            read = self.reads[i]
            part_s = float(read.attributes.get('S.ppm',0))/self.total_s
            part_e = float(read.attributes.get('E.ppm',0))/self.total_e
            ordered_ends = sorted([part_s, part_e])
            ratio = (ordered_ends[1]+0.01) / (ordered_ends[0]+0.01)
            average_ratio += ratio

        average_ratio = average_ratio / len(traceback)
        return average_ratio
    
    cpdef tuple get_name_for(self, index):
        """Given its relationship to the annotation set,
        generates a gene_id and transcript_id for
        row [index] of the membership matrix."""
        gene_id = 'bookend'
        transcript_id = 'bookend.1'
        return gene_id, transcript_id
    
    cpdef transcript_from_membership(self, int index):
        """Produces a collection of RNAseqMapping objects from the items in membership"""
        gene_id, transcript_id = self.get_name_for(index)
        e_tag = s_tag = True
        members = list(np.where(self.membership[index,:-4]==1)[0])
        ranges = []
        strand = self.strand_array[index]
        chrom = self.chrom
        source = 0
        
        elementAttributes = {}
        elementAttributes['gene_id'] = gene_id
        elementAttributes['transcript_id'] = transcript_id
        elementAttributes['source'] = 'bookend'
        elementAttributes['assemblies'] = self.rep_array[index]
        elementAttributes['assembled_in'] = ''
        weight = 0
        elementAttributes['S.reads'] = 0
        elementAttributes['S.capped'] = 0
        elementAttributes['E.reads'] = 0
        for i in self.traceback[index]:
            read = self.reads[i]
            weight += float(read.attributes['TPM'])
            elementAttributes['S.reads'] += float(read.attributes.get('S.reads', 0))
            elementAttributes['S.capped'] += float(read.attributes.get('S.capped', 0))
            elementAttributes['E.reads'] += float(read.attributes.get('E.reads', 0))
            elementAttributes['assembled_in'] += ',{}'.format(read.attributes.get('source',''))

        elementAttributes['assembled_in'] = elementAttributes['assembled_in'].lstrip(',')
        left = right = -1
        current_frag = (left, right)
        number_of_members = len(members)
        last = number_of_members - 1
        for i in range(number_of_members):
            left, right = self.frags[members[i]]
            if i == 0: # First member
                lbp = self.get_border_branchpoint(left, 0, strand)
                left = lbp.pos
                if strand == 1:
                    capped = lbp.percent_capped() >= self.cap_percent
            
            if i == last: # Last member
                rbp = self.get_border_branchpoint(right, 1, strand)
                right = rbp.pos
                if strand == -1:
                    capped = rbp.percent_capped() >= self.cap_percent
            
            if current_frag[1] == -1: # Uninitialized
                current_frag = (left, right)
            elif current_frag[1] == left: # Contiguous with the next member
                current_frag = (current_frag[0], right)
            else: # Jumped
                ranges.append(current_frag)
                current_frag = (left, right)
        
        ranges.append(current_frag)
        splice = [True] * (len(ranges)-1)
        elementData = ru.ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)
        readObject = ru.RNAseqMapping(elementData, elementAttributes)
        return readObject
    
    cpdef int matching_ref(self, read):
        """Returns the index number in ref_reads of the best-matching ref transcript."""
        cdef int match, smallest_diff, diff
        match = smallest_diff = -1
        for i in range(len(self.ref_reads)):
            ref = self.ref_reads[i]
            if read.overlaps(ref):
                if read.is_compatible(ref, ignore_ends=True, ignore_source=True):
                    diff = read.diff(ref)
                    if smallest_diff == -1 or diff < smallest_diff:
                        smallest_diff = diff
                        match = i
        
        return match
    
    cpdef int sense_ref(self, read):
        """Returns the index number in ref_reads of the best-matching ref transcript."""
        cdef int match, smallest_diff, diff
        match = smallest_diff = -1
        for i in range(len(self.ref_reads)):
            ref = self.ref_reads[i]
            if read.sense_match(ref, 1):
                diff = read.diff(ref)
                if smallest_diff == -1 or diff < smallest_diff:
                    smallest_diff = diff
                    match = i
        
        return match

    cpdef int antisense_ref(self, read):
        """Returns the index number in ref_reads of the best-matching ref transcript."""
        cdef int match, smallest_diff, diff
        match = smallest_diff = -1
        for i in range(len(self.ref_reads)):
            ref = self.ref_reads[i]
            if read.antisense_match(ref, self.end_extend):
                diff = read.diff(ref)
                if smallest_diff == -1 or diff < smallest_diff:
                    smallest_diff = diff
                    match = i
        
        return match
    

##########################################


cpdef np.ndarray get_information_content(np.ndarray[char, ndim=2] membership_matrix):
    """Given a matrix of membership values,
    return sum(abs(frag_membership)) for each row.
    Maximum information is number_of_frags + 4"""
    return np.sum(np.abs(membership_matrix), 1, dtype=np.int32)

cpdef np.ndarray get_member_content(np.ndarray[char, ndim=2] membership_matrix):
    """Given a matrix of membership values,
    return sum(abs(frag_membership)) for each row.
    Maximum information is number_of_frags + 4"""
    return np.sum(membership_matrix==1, 1, dtype=np.int32)

cpdef np.ndarray calculate_overlap(np.ndarray[char, ndim=2] membership_matrix, np.ndarray information_content, np.ndarray strand_array):
    """Given a matrix of membership values (1, 0, or -1; see locus.build_membership_matrix),
    output a new read x read square matrix with a overlap code.
    For each (a,b) pair:
        -1 - a is incompatible with b (at least one (-1,1) pair)
        0 - no overlapping memberships between a and b (no (1,1) pairs or (-1,1) pairs)
        1 - a is partially matched to b (at least on (1,1) pair, no (-1,1) pairs)
        2 - a is contained within b (all 1's of a are 1 in b)
    
    Compatibility is not commutative: comp(a,b) is not necessarily == comp(b,a).
    If an information_content vertex is not supplied, it will be calculated by the
    sum(abs(frag_membership)) of each row.
    """
    cdef int ia, ib, shared
    cdef char sa, sb, s
    cdef bint overlapping
    cdef Py_ssize_t a, b, i, maxlen, number_of_reads, number_of_frags
    
    number_of_reads = membership_matrix.shape[0]
    number_of_frags = membership_matrix.shape[1]
    if len(information_content) != number_of_reads:
        information_content = get_information_content(membership_matrix)
    
    cdef np.ndarray[char, ndim=2] overlap_matrix = np.zeros((number_of_reads, number_of_reads), dtype=np.int8) # Container for overlap information
    cdef char [:] STRAND_ARRAY = strand_array
    cdef int [:] INFO = information_content
    cdef char [:, :] MEMBERSHIP = membership_matrix
    cdef char [:, :] COMPATIBILITY = overlap_matrix
    maxlen = COMPATIBILITY.shape[0]

    for a in range(maxlen):
        for b in range(a,maxlen):
            if b == a:
                COMPATIBILITY[a,b] = 2 # A read is necessarily identical to itself
                continue
            
            sa = STRAND_ARRAY[a]
            sb = STRAND_ARRAY[b]
            if sa != 0 and sb != 0:
                if sa != sb: # Reads are on different strands, 
                    COMPATIBILITY[a,b] = -1
                    COMPATIBILITY[b,a] = -1
                    continue
            
            if sa > 0 or sb > 0:
                s = 1
            elif sa < 0 or sb < 0:
                s = -1
            else:
                s = 0
            
            shared = 0
            overlapping = False
            for i in range(MEMBERSHIP.shape[1]):
                ia = MEMBERSHIP[a,i]
                ib = MEMBERSHIP[b,i]
                if ia == 0 or ib == 0: # Incomplete information
                    continue
                
                if ia != ib: # A (-1, 1) pair exists
                    COMPATIBILITY[a,b] = -1
                    COMPATIBILITY[b,a] = -1
                    shared = -1
                    break # No need to evaluate the rest of the frags
                
                if overlapping == False:
                    if ia == 1 and ib == 1:
                        overlapping = True
                
                shared += 1
                
            if shared <= 0:
                continue
            else:
                if shared == INFO[a]: # All information of a is shared in b
                    COMPATIBILITY[a,b] = 2
                elif overlapping and s >= 0:
                    COMPATIBILITY[a,b] = 1
                
                if shared == INFO[b]: # All information of b is shared in a
                    COMPATIBILITY[b,a] = 2
                elif overlapping and s <= 0:
                    COMPATIBILITY[b,a] = 1
    
    return overlap_matrix

cpdef bint passes_threshold(np.ndarray array, int max_gap, float threshold=1):
    """Returns boolean of whether a contiguous region of values
    less than threshold and longer than max_gap exists in array."""
    cdef Py_ssize_t i
    cdef int run = 0
    cdef float a
    cdef bint passed = False
    for i in range(array.shape[0]):
        a = array[i]
        if a < threshold:
            run += 1
            if run > max_gap:
                return False
        else:
            run = 0
            if passed is False:
                passed = True
    
    return passed

cpdef np.ndarray lerp(np.ndarray array, float null_value):
    """Given an array with some null values, in-fill null ranges by
    linear interpolation from the flanking values."""


cpdef np.ndarray remove_ends(np.ndarray[char, ndim=2] membership_matrix):
    """Given a full membership matrix, return a reduced version in which
    end information is ignored."""
    cdef:
        np.ndarray endless_matrix, boundaries
        Py_ssize_t i, columns
        int first, last
    
    endless_matrix = copy.copy(membership_matrix[:,:-4])
    boundaries = np.apply_along_axis(first_and_last, 1, endless_matrix)
    columns = endless_matrix.shape[1]
    for i in range(endless_matrix.shape[0]):
        first, last = boundaries[i,:2]
        if first >= 0:
            endless_matrix[i,:first] = 0
        
        if last >= 0:
            endless_matrix[i,(last+1):] = 0
    
    return endless_matrix

cpdef (int, int) first_and_last(np.ndarray[char, ndim=1] membership_row):
    """Returns a (left,right) tuple of the first and last member positions."""
    cdef np.ndarray indices = np.where(membership_row==1)[0]
    if len(indices) > 0:
        return (indices[0],indices[-1])
    else:
        return (-1, -1)

cpdef list find_breaks(np.ndarray[char, ndim=2] membership_matrix, bint ignore_ends=True):
    """Identifies all points along the membership array where it could
    be cleanly divided in two, with reads entirely on one side or the other."""
    cdef np.ndarray boundaries, gaps, boolarray
    cdef list breaks = []
    if ignore_ends:
        boundaries = np.apply_along_axis(first_and_last, 1, membership_matrix[:,:-4])
        gaps = np.where(np.sum(membership_matrix[:,:-4]==1,0)==0)[0]
    else:
        boundaries = np.apply_along_axis(first_and_last, 1, membership_matrix)
        gaps = np.where(np.sum(membership_matrix==1,0)==0)[0]
    
    for g in gaps:
        boolarray = g > boundaries
        if not np.any(np.sum(boolarray, 1)==1): # G isn't inside any membership ranges
            if len(np.unique(boolarray)) == 2:
                breaks.append(g)
    
    return breaks

cpdef np.ndarray find_linear_chains(np.ndarray[char, ndim=2] overlap_matrix, np.ndarray resolve_order):
    cdef:
        np.ndarray edges, ingroup, outgroup, in_chain, putative_chain_starts, containment
        int chain, v, next_v, c
        set contained, parent_chains
    
    edges = overlap_matrix == 1
    edges = np.logical_and(edges, overlap_matrix.transpose()!=2)
    ingroup = np.sum(edges, axis=0, dtype=np.int32)
    outgroup = np.sum(edges, axis=1, dtype=np.int32)
    containment = overlap_matrix == 2
    np.put(containment, range(0,containment.shape[0]**2,containment.shape[0]+1), False, mode='wrap') # Blank out the diagonal (self-containments)
    contained = set(np.where(np.sum(containment,axis=1)>0)[0])
    in_chain = np.zeros(len(ingroup), dtype=np.int32)
    putative_chain_starts = np.where(outgroup == 1)[0]
    chain = 0
    for v in putative_chain_starts:
        if in_chain[v] == 0: # v is unvisited
            chain += 1
            next_v  = np.where(edges[v,:])[0][0]
            if ingroup[next_v] == 1:
                in_chain[v] = chain
            
            while ingroup[next_v] == 1:
                if in_chain[next_v] != 0: # next_v is already visited, join v to this chain
                    in_chain[v] = in_chain[next_v]
                    break
                
                in_chain[next_v] = chain
                if outgroup[next_v] != 1:break
                next_v = np.where(edges[next_v,:])[0][0]
                if in_chain[next_v] == chain: break
    
    # If all a contained element's containers are in the same chain, add it to this chain
    for c in resolve_order:
        if c in contained:
            parent_chains = set(in_chain[containment[c,:]])
            if len(parent_chains) == 1:
                in_chain[c] = parent_chains.pop()
    
    return in_chain


def sum_subset(mask, array_to_mask):
    """Given an array and a boolean mask, return
    the sum of array values at which mask was True"""
    return array_to_mask[mask].sum()

