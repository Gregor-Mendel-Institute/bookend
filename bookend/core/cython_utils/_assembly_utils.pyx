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
    cdef public float weight, capped
    cdef public str tag
    cdef public object positions
    """Represents a reference point for a Start or End site."""
    def __init__(self, left, right, peak, weight, endtype):
        self.left, self.right, self.peak, self.weight, self.endtype = left, right, peak, weight, endtype
        self.tag, self.strand = [('S', 1), ('E', 1), ('S', -1), ('E', -1)][self.endtype]
        self.positions = Counter()
        if self.endtype in [0, 3]:
            self.terminal = left
        else:
            self.terminal = right
    
    def __repr__(self):
        strand = ['.','+','-'][self.strand]
        return '{}{}{} ({})'.format(self.tag, strand, self.peak, self.weight)
    
    def __eq__(self, other): return self.span() == other.span()
    def __ne__(self, other): return self.span() != other.span()
    def __gt__(self, other): return self.span() >  other.span()
    def __ge__(self, other): return self.span() >= other.span()
    def __lt__(self, other): return self.span() <  other.span()
    def __le__(self, other): return self.span() <= other.span()
    
    def add(self, int position, float weight, bint capped=False):
        self.positions[position] += weight
        self.capped += weight*capped
    
    cpdef (int, int) span(self):
        return (min(self.positions.keys()), max(self.positions.keys()))
    
    def most_common(self):
        return self.positions.most_common(1)[0][0]
    
    def write_as_bed(self, chrom):
        print('{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, self.left, self.right, self.tag, self.weight, {-1:'-',1:'+',0:'.'}[self.strand]))

cdef class Locus:
    cdef public int chrom, leftmost, rightmost, extend, end_extend, number_of_elements, min_overhang, chunk_number, oligo_len, min_intron_length
    cdef public bint naive, allow_incomplete, use_attributes, ignore_ends, require_cap, splittable, verbose
    cdef public tuple reads, frags
    cdef public float weight, bases, raw_bases, minimum_proportion, cap_bonus, cap_filter, intron_filter, antisense_filter, dead_end_penalty
    cdef public dict J_plus, J_minus, end_ranges, source_lookup, adj, exc
    cdef public set branchpoints, SPbp, EPbp, SMbp, EMbp
    cdef public list transcripts, traceback, sources, subproblem_indices, splits
    cdef public object graph
    cdef EndRange nullRange
    cdef Locus sublocus
    cdef public np.ndarray depth_matrix, strandratio, cov_plus, cov_minus, depth, read_lengths, discard_frags, member_lengths, frag_len, frag_by_pos, strand_array, weight_array, rep_array, membership, overlap, information_content, member_content, frag_strand_ratios, member_weights
    def __init__(self, chrom, chunk_number, list_of_reads, max_gap=50, end_cluster=200, min_overhang=3, reduce=True, minimum_proportion=0.01, min_intron_length=50, antisense_filter=0.01, cap_bonus=5, cap_filter=.1, complete=False, verbose=False, naive=False, intron_filter=0.10, use_attributes=False, oligo_len=20, ignore_ends=False, allow_incomplete=False, require_cap=False, splittable=True):
        self.nullRange = EndRange(-1, -1, -1, -1, -1)
        self.oligo_len = oligo_len
        self.transcripts = []
        self.traceback = []
        self.branchpoints = set()
        self.chunk_number = chunk_number
        self.naive = naive
        self.minimum_proportion = minimum_proportion
        self.min_intron_length = min_intron_length
        self.intron_filter = intron_filter
        self.antisense_filter = antisense_filter
        self.min_overhang = min_overhang
        self.chrom = chrom
        self.cap_bonus = cap_bonus
        self.cap_filter = cap_filter
        self.use_attributes = use_attributes
        self.allow_incomplete = allow_incomplete
        self.ignore_ends = ignore_ends
        self.require_cap = require_cap
        self.splittable = splittable
        self.verbose = verbose
        if self.ignore_ends:
            self.dead_end_penalty = 1
        elif self.allow_incomplete:
            self.dead_end_penalty = 0.1
        else:
            self.dead_end_penalty = 0
        
        if len(list_of_reads) > 0:
            self.leftmost, self.rightmost = ru.range_of_reads(list_of_reads)
            if self.ignore_ends:
                self.allow_incomplete = True
                for read in list_of_reads:
                    read.s_tag, read.e_tag = False, False
                    if len(read.splice)==0:read.strand = 0
            else: # Don't bother processing a list of reads without at least one same-stranded end pair
                if not ru.has_ends(list_of_reads, self.require_cap):
                    return
            
            self.reads = tuple(list_of_reads) # Cannot be mutated            
            self.read_lengths = np.array([r.get_length()+self.oligo_len*r.s_tag+self.oligo_len*r.e_tag for r in self.reads], dtype=np.int32)
            self.raw_bases = np.sum(self.read_lengths * np.array([read.weight for read in self.reads]))
            self.sources = ru.get_sources(list_of_reads)
            self.source_lookup = ru.get_source_dict(self.sources)
            self.extend = max_gap
            self.end_extend = end_cluster
            self.depth_matrix, self.J_plus, self.J_minus = ru.build_depth_matrix(self.leftmost, self.rightmost, self.reads, self.use_attributes)
            Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
            covstranded = np.sum(self.depth_matrix[(covp,covm),:],axis=0)
            strandedpositions = np.where(covstranded > 0)[0]
            if strandedpositions.shape[0] > 0: # Some reads to inform the strand
                strandratio = np.array(np.interp(range(covstranded.shape[0]), strandedpositions, self.depth_matrix[covp,strandedpositions]/covstranded[strandedpositions]),dtype=np.float32)
                strandratio[strandratio < self.antisense_filter] = 0
                strandratio[strandratio > 1-self.antisense_filter] = 1
                self.strandratio = strandratio
                self.cov_plus = self.depth_matrix[covp,:] + self.depth_matrix[covn,:]*self.strandratio
                self.cov_minus = self.depth_matrix[covm,:] + self.depth_matrix[covn,:]*(1-self.strandratio)
                self.depth = self.cov_plus + self.cov_minus
            else:
                self.strandratio = np.full(self.depth_matrix.shape[1], .5)
                self.cov_plus = self.depth_matrix[covn,:]*self.strandratio
                self.cov_minus = self.depth_matrix[covn,:]*self.strandratio
                self.depth = self.depth_matrix[covn,:]
            
            self.prune_junctions(self.min_intron_length)
            self.generate_branchpoints()
            # Split locus into coherent subchunks
            self.splits = []
            if self.splittable:
                self.splits = self.split_chunk()
            
            if len(self.splits) > 0:
                if self.verbose:print('Processing {} subchunks'.format(len(self.splits)+1))
                for subchunk in ru.generate_subchunks(list_of_reads, self.splits):
                    self.chunk_number += 1
                    sublocus = Locus(self.chrom, self.chunk_number, subchunk, max_gap=self.extend, end_cluster=self.end_extend, min_overhang=self.min_overhang, reduce=True, minimum_proportion=self.minimum_proportion, min_intron_length=self.min_intron_length, antisense_filter=self.antisense_filter, cap_bonus=self.cap_bonus, cap_filter=self.cap_filter, complete=False, verbose=False, naive=self.naive, intron_filter=self.intron_filter, use_attributes=self.use_attributes, oligo_len=self.oligo_len, ignore_ends=self.ignore_ends, allow_incomplete=self.allow_incomplete, require_cap=self.require_cap, splittable=False)
                    self.transcripts += sublocus.transcripts
            else:
                self.build_membership_matrix()
                if self.membership.shape[0] > 0:
                    self.build_overlap_matrix()
                    self.build_graph(reduce)
                
                if self.bases > 0:
                    self.assemble_transcripts()
    
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
    
    cpdef list split_chunk(self):
        """Check whether, given the filtered splice junctions, the list_of_reads can be separated
        into a list of smaller coherent read sets with no overlap.
        Returns the list of read indices that begin new coherent subchunks."""
        cdef np.ndarray split_bool, difs, run_starts, run_ends
        cdef str k, kl, kr
        cdef int rs, rl
        split_bool = self.depth == 0
        for k in self.J_plus.keys():
            kl,kr = k.split(':')
            split_bool[int(kl)-1:int(kr)+1] = False
        
        for k in self.J_minus.keys():
            kl,kr = k.split(':')
            split_bool[int(kl)-1:int(kr)+1] = False
        
        difs = np.diff(np.hstack(([0], split_bool, [0])))
        run_starts = np.where(difs > 0)[0]
        run_ends = np.where(difs < 0)[0]
        return [rs+self.leftmost for rs,rl in zip(run_starts, run_ends-run_starts) if rl > self.extend]
    
    cpdef void prune_junctions(self, int min_intron_length=50):
        """If a splice junction represents < minimum_proportion of junction-spanning reads
        that either share its donor or its acceptor site, remove it."""
        cdef dict jdict, leftsides, rightsides
        cdef list keys, spans
        cdef str k
        cdef int i,l,r, Sp, Ep, Sm, Em, covp, covm, covn
        cdef set prune
        cdef np.ndarray counts
        cdef float total, jcov, spanning_cov
        cdef bint passes
        Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
        for jdict in [self.J_plus, self.J_minus]:
            if jdict:
                prune = set()
                keys = sorted(list(jdict.keys()))
                spans = [self.string_to_span(k) for k in keys]
                leftsides = {}
                rightsides = {}
                for i in range(len(keys)):
                    l,r = spans[i]
                    spanning_cov = np.max(self.depth[l:r+1])
                    jcov = jdict[keys[i]]
                    if jcov < spanning_cov * self.minimum_proportion or r-l < min_intron_length:
                        prune.add(keys[i])
                    else:
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
            float threshold_depth, cumulative_depth, cutoff
            int l, r, p, Sp, Ep, Sm, Em, covp, covm, covn
            (int, int) span
            EndRange rng
            set prohibited_plus, prohibited_minus, prohibited_positions, bpset
            list gaps
        
        Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
        self.branchpoints = set()
        self.SPbp, self.EPbp, self.SMbp, self.EMbp = set(), set(), set(), set()
        self.end_ranges = dict()
        prohibited_plus, prohibited_minus = set(), set()
        for j in self.J_plus.keys():
            span = self.string_to_span(j)
            self.branchpoints.update(list(span))
            if self.min_overhang > 0:
                prohibited_plus.update(range(span[0]-self.min_overhang, span[0]+self.min_overhang+1))
                prohibited_plus.update(range(span[1]-self.min_overhang, span[1]+self.min_overhang+1))
        
        for j in self.J_minus.keys():
            span = self.string_to_span(j)
            self.branchpoints.update(list(self.string_to_span(j)))
            if self.min_overhang > 0:
                prohibited_minus.update(range(span[0]-self.min_overhang, span[0]+self.min_overhang+1))
                prohibited_minus.update(range(span[1]-self.min_overhang, span[1]+self.min_overhang+1))
        
        for endtype in [Ep, Em, Sp, Sm]:
            if endtype == Sp:
                bpset = self.SPbp
                pos = np.where(np.sum(self.depth_matrix[[Sp,Cp],],axis=0)>0)[0]
                vals = np.power(self.depth_matrix[Sp, pos] + self.cap_bonus*self.depth_matrix[Cp, pos],2)/self.cov_plus[pos]*(self.strandratio[pos]!=0).astype(float)
                prohibited_positions = prohibited_plus
            elif endtype == Sm:
                bpset = self.SMbp
                pos = np.where(np.sum(self.depth_matrix[[Sm,Cm],],axis=0)>0)[0]
                vals = np.power(self.depth_matrix[Sm, pos] + self.cap_bonus*self.depth_matrix[Cm, pos],2)/self.cov_minus[pos]*(self.strandratio[pos]!=1).astype(float)
                prohibited_positions = prohibited_minus
            elif endtype == Ep:
                bpset = self.EPbp
                pos = np.where(self.depth_matrix[endtype,]>0)[0]
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_plus[pos]*(self.strandratio[pos]!=0).astype(float)
                prohibited_positions = prohibited_plus
            elif endtype == Em:
                bpset = self.EMbp
                pos = np.where(self.depth_matrix[endtype,]>0)[0]
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_minus[pos]*(self.strandratio[pos]!=1).astype(float)
                prohibited_positions = prohibited_minus
            
            self.end_ranges[endtype] = self.make_end_ranges(pos, vals, endtype, prohibited_positions)
            if endtype in [Sp, Sm]:
                self.resolve_overlapping_ends(endtype)
                self.enforce_cap_filter(endtype)
        
        for endtype in [Ep, Em, Sp, Sm]:
            for rng in self.end_ranges[endtype]:
                self.branchpoints.add(rng.terminal)
                bpset.add(rng.terminal)
                
        self.branchpoints.add(0)
        self.branchpoints.add(len(self))
    
    cpdef void enforce_cap_filter(self, int endtype):
        cdef int caps
        cdef np.ndarray allstarts
        if not self.require_cap:
            return
        
        if endtype == 0:
            caps = 4
        elif endtype == 2:
            caps = 5
        else:
            return
        
        allstarts = self.depth_matrix[endtype,:]+self.depth_matrix[caps,:]
        self.end_ranges[endtype] = [
            rng for rng in self.end_ranges[endtype]
            if np.sum(self.depth_matrix[caps, rng.left:rng.right]) >= np.sum(allstarts[rng.left:rng.right])*self.cap_filter
        ]
        return
    
    cpdef void resolve_overlapping_ends(self, int endtype):
        """If same-stranded start and end clusters overlap, split
        to avoid malformed entries where start is downstream of end.
        """
        cdef EndRange LR, RR
        cdef int ltype, rtype
        cdef list newlefts, newrights
        cdef set added_lefts, added_rights
        newlefts, newrights = [], []
        added_lefts, added_rights = set(), set()
        if endtype == 0:
            ltype = 0
            rtype = 1
            lefts = self.end_ranges[ltype]
            rights = self.end_ranges[rtype]
        elif endtype == 2:
            ltype = 3
            rtype = 2
            lefts = self.end_ranges[ltype]
            rights = self.end_ranges[rtype]
        
        for LR in lefts:
            for RR in rights:
                if RR.peak not in added_rights: # 
                    if RR.right < LR.left:
                        newrights += [RR]
                        added_rights.add(RR.peak)
                        continue
                    elif RR.left > LR.right:
                        newlefts += [LR]
                        added_lefts.add(LR.peak)
                        continue
                    
                    if LR.right > RR.left and LR.left < RR.right: # overlapping
                        if LR.peak > RR.peak: # out-of-order peaks
                            left_positions = Counter({k:v for k,v in LR.positions.items() if k < RR.peak})
                            if len(left_positions) > 0:
                                splitRangeLeft = EndRange(min(left_positions), max(left_positions), left_positions.most_common(1)[0][0], sum(left_positions.values()), ltype)
                                splitRangeLeft.positions = left_positions
                                newlefts += [splitRangeLeft]
                            
                            right_positions = Counter({k:v for k,v in LR.positions.items() if k > RR.peak})
                            if len(right_positions) > 0:
                                splitRangeRight = EndRange(min(right_positions), max(right_positions), right_positions.most_common(1)[0][0], sum(right_positions.values()), ltype)
                                splitRangeRight.positions = right_positions
                                newlefts += [splitRangeRight]
                                if RR.right > splitRangeRight.peak: # Also split RR
                                    left_positions = Counter({k:v for k,v in RR.positions.items() if k < splitRangeRight.peak})
                                    right_positions = Counter({k:v for k,v in RR.positions.items() if k > splitRangeRight.peak})
                                    if len(left_positions) > 0:
                                        splitRangeLeft = EndRange(min(left_positions), max(left_positions), left_positions.most_common(1)[0][0], sum(left_positions.values()), rtype)
                                        splitRangeLeft.positions = left_positions
                                        newrights += [splitRangeLeft]
                                    
                                    if len(right_positions) > 0:
                                        splitRangeRight = EndRange(min(right_positions), max(right_positions), right_positions.most_common(1)[0][0], sum(right_positions.values()), rtype)
                                        splitRangeRight.positions = right_positions
                                        newrights += [splitRangeRight]
                                    
                                else:
                                    newrights += [RR]
                                    added_rights.add(RR.peak)
                            
                            added_lefts.add(LR.peak)
                            added_rights.add(RR.peak)
                        else: # Peaks are in the correct order, do nothing
                            newlefts += [LR]
                            added_lefts.add(LR.peak)
            
            if LR.peak not in added_lefts:
                newlefts += [LR]
                added_lefts.add(LR.peak)
        
        for RR in rights:
            if RR.peak not in added_rights:
                newrights += [RR]
                added_rights.add(RR.peak)
        
        self.end_ranges[ltype] = newlefts
        self.end_ranges[rtype] = newrights
    
    cpdef bint passes_cap_filter(self, EndRange rng):
        """"Checks if an EndRange should be kept. Sp/Sm EndRanges with >50% flowthrough
        must be >self.cap_filter uuG to be retained
        """
        cdef float cov_in, cov_through, caps, noncaps
        cdef int Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn
        cdef int bp, closest_bp
        Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
        if rng.endtype == Sp:
            if rng.left == 0:
                return True
            
            closest_bp = max([0]+[bp for bp in self.branchpoints if bp < rng.right])
            if closest_bp == 0:
                return True
            
            cov_through = np.mean(self.cov_plus[closest_bp:rng.right])
            cov_in = self.cov_plus[rng.right]
            if cov_through < cov_in*self.intron_filter:
                return True
            
            noncaps = np.sum(self.depth_matrix[Sp, rng.left:rng.right])
            caps = np.sum(self.depth_matrix[Cp, rng.left:rng.right])
            if caps >= self.cap_filter * (caps + noncaps):
                return True
            
            return False
        elif rng.endtype == Sm:
            if rng.right >= self.cov_minus.shape[0]:
                return True
            
            closest_bp = min([self.cov_minus.shape[0]]+[bp for bp in self.branchpoints if bp > rng.left+1])
            if closest_bp == self.cov_minus.shape[0]:
                return True
            
            cov_through = np.mean(self.cov_minus[rng.left+1:closest_bp])
            cov_in = self.cov_minus[rng.left]
            if cov_through < cov_in*self.intron_filter:
                return True
            
            noncaps = np.sum(self.depth_matrix[Sm, rng.left:rng.right])
            caps = np.sum(self.depth_matrix[Cm, rng.left:rng.right])
            if caps >= self.cap_filter * (caps + noncaps):
                return True
            
            return False
        else:
            return True
    
    cpdef list make_end_ranges(self, np.ndarray pos, np.ndarray vals, int endtype, set prohibited_positions):
        """Returns a list of tuples that (1) filters low-signal positions
        and (2) clusters high-signal positions within self.end_extend.
        Returns list of (l,r) tuples demarking the edges of clustered end positions."""
        cdef:
            np.ndarray value_order, passes_threshold
            EndRange e
            int p, maxp, prohibit_pos, i
            float cumulative, threshold, v, maxv, weight
            list filtered_pos, end_ranges
            (int, int) current_range
            bint passed_prohibit
            dict positions
        
        if len(pos) == 0:
            return []
        elif len(pos) == 1:
            return [EndRange(pos[0], pos[0]+1, pos[0], vals[0], endtype)]
        
        passes_threshold = vals > self.minimum_proportion*np.sum(vals)
        if endtype in [0, 3]:
            passes_threshold[0] = True
        elif endtype in [1, 2]:
            passes_threshold[-1] = True
        
        if np.sum(passes_threshold) == 0:
            return []
        
        filtered_pos = [(p,v) for p,v in zip(pos[passes_threshold], vals[passes_threshold]) if p not in prohibited_positions]
        if len(filtered_pos) == 0:
            return []
        
        p,v = filtered_pos[0]
        maxv = v
        maxp = p
        weight = v
        positions = {p:v}
        current_range = (p, p+1)
        end_ranges = []
        prohibited = iter(sorted(prohibited_positions))
        passed_prohibit = False
        prohibit_pos = -1
        while prohibit_pos < p: # Initialize prohibit_pos as the first prohibited_position right of current_range 
            try:
                prohibit_pos = next(prohibited)
            except StopIteration:
                prohibit_pos = -1
                break
        
        for p,v in filtered_pos[1:]: # Iterate over all positions that passed the threshold
            passed_prohibit = prohibit_pos > -1 and p > prohibit_pos
            if passed_prohibit or p - self.end_extend > current_range[1]:  # Must start a new range
                e = EndRange(current_range[0], current_range[1], maxp, weight, endtype)
                e.positions = Counter(positions)
                if self.passes_cap_filter(e):
                    end_ranges.append(e)
                
                current_range = (p, p+1)
                positions = {p:v}
                maxp = p
                maxv = v
                weight = v
                if passed_prohibit: # Push the prohibit_pos past current p
                    passed_prohibit = False
                    while prohibit_pos < p:
                        try:
                            prohibit_pos = next(prohibited)
                        except StopIteration:
                            prohibit_pos = -1
                            break
            else:# Can continue the last range
                current_range = (current_range[0], p+1)
                positions[p] = v
                weight += v
                if v > maxv or (v == maxv and endtype in [1,2]):
                    maxv, maxp = v, p
        
        e = EndRange(current_range[0], current_range[1], maxp, weight, endtype)
        e.positions = Counter(positions)
        if self.passes_cap_filter(e):
            end_ranges.append(e)
        
        return end_ranges
    
    cpdef int end_of_cluster(self, int pos, float weight, list end_ranges, int extend, bint capped):
        """Returns the terminal position of the EndRange object that
        contains pos, if one exists. Else returns -1"""
        cdef EndRange rng
        rng = self.get_end_cluster(pos, weight, end_ranges, extend, capped)
        return rng.terminal
    
    cpdef EndRange get_end_cluster(self, int pos, float weight, list end_ranges, int extend, bint capped=False):
        """Returns the most common position of the EndRange object that
        contains pos, if one exists. Else returns -1"""
        cdef EndRange rng, bestrng
        cdef int dist, bestdist
        bestdist = -1
        bestrng = self.nullRange
        for rng in end_ranges:
            if pos >= max(0, rng.left-extend) and pos <= min(self.frag_by_pos.shape[0], rng.right+extend):
                dist = abs(pos-rng.peak)
                if bestdist == -1 or dist < bestdist:
                    bestrng = rng
                    bestdist = dist
        
        bestrng.add(pos, weight, capped)
        return bestrng
    
    cpdef str span_to_string(self, (int, int) span):
        """Converts a tuple of two ints to a string connected by ':'"""
        return '{}:{}'.format(span[0], span[1])
    
    cpdef (int, int) string_to_span(self, str string):
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
        values are produced for each nucleotide of the object as:
            self.frag_by_pos - the unique fragment that each nucleotide belongs to
        After this, the table is compressed to produce:
            self.reduced_membership - a table of unique membership patterns
            self.membership_weights - (1 per row of reduced_membership) number of reads
        """
        cdef Py_ssize_t a, b, i, number_of_reads, number_of_frags
        cdef int membership_width, read_length, lost_ends, number_of_members
        cdef np.ndarray[char, ndim=1] membership
        cdef list temp_frags, bp_positions, hashes
        cdef (int, int) block, span
        cdef str membership_hash, null_hash
        cdef dict reps, source_weights, member_weights, membership_lengths, unhash
        unhash = {'_':-1, ' ':0, '*':1}
        reps, source_weights, member_weights, membership_lengths = {}, {}, {}, {}
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
            else:
                self.frag_strand_ratios[i] = .5 # Naive if no known stranded coverage
        
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
        membership_width = number_of_frags+4
        null_hash = '_'*membership_width
        self.discard_frags = self.apply_intron_filter(threshold)
        for i in range(number_of_reads): # Read through the reads once, cataloging frags and branchpoints present/absent
            read = self.reads[i]
            membership = self.calculate_membership(membership_width, read.ranges, read.splice, read.strand, read.s_tag, read.e_tag, read.capped, read.weight)
            membership_hash = ''.join([[' ', '*', '_'][m] for m in membership])
            if membership_hash != null_hash: # Read wasn't discarded
                if read.strand == 1:
                    lost_ends = sum([read.s_tag and membership[membership_width-4]!=1, read.e_tag and membership[membership_width-3]!=1])
                elif read.strand == -1:
                    lost_ends = sum([read.s_tag and membership[membership_width-2]!=1, read.e_tag and membership[membership_width-1]!=1])
                else:
                    lost_ends = 0
                
                read_length = self.read_lengths[i] - self.oligo_len*lost_ends
                if membership_hash in reps.keys(): # This pattern has already been seen
                    reps[membership_hash] += read.weight
                    source_weights[membership_hash][self.source_lookup[read.source]] += read.weight * read_length / membership_lengths[membership_hash]
                else: # This is the first time encountering this pattern
                    reps[membership_hash] = read.weight
                    membership_lengths[membership_hash] = np.sum(self.frag_len[membership == 1])
                    source_weights[membership_hash] = np.zeros(len(self.source_lookup), dtype=np.float32)
                    source_weights[membership_hash][self.source_lookup[read.source]] += read.weight * read_length / membership_lengths[membership_hash]
        
        hashes = sorted(reps.keys())
        number_of_members = len(hashes)
        self.member_lengths = np.zeros(number_of_members, dtype=np.int32)
        self.rep_array = np.zeros(number_of_members, dtype=np.float32)
        self.membership = np.zeros((number_of_members, membership_width), dtype=np.int8)
        self.member_weights = np.zeros((number_of_members,membership_width), dtype=np.float32)
        self.weight_array = np.zeros((number_of_members,len(self.source_lookup)), dtype=np.float32)
        self.strand_array = np.zeros(number_of_members, dtype=np.int8)
        for i in range(number_of_members):
            membership_hash = hashes[i]
            self.membership[i,:] = [unhash[h] for h in membership_hash]
            self.weight_array[i,:] = source_weights[membership_hash]
            self.rep_array[i] = reps[membership_hash]
            self.member_weights[i,self.membership[i,:]==1] = np.sum(self.weight_array[i,:])
            self.member_weights[i,self.membership[i,:]==-1] = self.rep_array[i]
            self.strand_array[i] = 0 - int(membership_hash[-4:-2]=='__') + int(membership_hash[-2:]=='__')
        
        if self.naive:
            self.weight_array = np.sum(self.weight_array,axis=1,keepdims=True)
        
        self.reduce_membership()
        self.filter_by_reps(threshold)
        self.filter_members_by_strand()
        self.reduce_membership()
        self.weight = np.sum(self.weight_array)
        self.number_of_elements = self.membership.shape[0]
        self.information_content = get_information_content(self.membership)
        self.member_content = get_member_content(self.membership)
        self.bases = np.sum(np.sum(self.weight_array, axis=1)*self.member_lengths)
        return False
    
    cpdef np.ndarray[char, ndim=1] calculate_membership(self, int width, list ranges, list splice, char strand, bint s_tag, bint e_tag, bint capped, float weight):
        """Given an RNAseqMapping object, defines a membership string that describes
        which frags are included (*), excluded (_) or outside of ( ) the read."""
        cdef np.ndarray[char, ndim=1] membership
        cdef np.ndarray members
        cdef Py_ssize_t j, source_plus, source_minus, sink_plus, sink_minus
        cdef (int, int) block, span
        cdef int last_rfrag, l, r, tl, tr, lfrag, rfrag
        cdef char s, jstrand
        cdef bint junctions_are_linear
        cdef list spans, splice_sites, sorted_splice_sites, intervening_junctions, skipped_frags
        cdef str junction_hash, junction, membership_hash
        membership = np.zeros(width, dtype=np.int8)
        source_plus, sink_plus, source_minus, sink_minus = range(width-4, width)
        last_rfrag = 0
        if strand == 1:
            membership[source_minus] = -1 # Read cannot have minus-stranded features
            membership[sink_minus] = -1 # Read cannot have minus-stranded features
        elif strand == -1:
            membership[source_plus] = -1 # Read cannot have plus-stranded features
            membership[sink_plus] = -1 # Read cannot have plus-stranded features
        
        for j in range(len(ranges)): # Run through each block range of read
            block = ranges[j]
            l = block[0] - self.leftmost
            r = block[1] - self.leftmost
            # Get membership information about the block
            lfrag = self.frag_by_pos[l]
            rfrag = self.frag_by_pos[r-1]
            if self.min_overhang > 0:
                if self.frag_len[lfrag] > self.min_overhang and l+self.min_overhang < len(self):
                    lfrag = self.frag_by_pos[l+self.min_overhang-1]
                
                if self.frag_len[rfrag] > self.min_overhang and r-self.min_overhang >= 0:
                    rfrag = self.frag_by_pos[r-self.min_overhang]
            
            if j == 0: # Starting block
                if strand == 1 and s_tag: # Left position is a 5' end
                    # Sp, Ep, Sm, Em
                    tl = self.end_of_cluster(l, weight, self.end_ranges[0], self.end_extend*capped, capped)
                    if tl >= 0: # 5' end is in an EndRange
                        membership[source_plus] = 1 # Add s+ to the membership table
                        l = tl
                        lfrag = self.frag_by_pos[l]
                        membership[0:lfrag] = -1 # Read cannot extend beyond source
                elif strand == -1 and e_tag: # Left position is a 3' end
                    tl = self.end_of_cluster(l, weight, self.end_ranges[3], self.end_extend, False)
                    if tl >= 0:
                        membership[sink_minus] = 1 # Add t- to the membership table
                        l = tl
                        lfrag = self.frag_by_pos[l]
                        membership[0:lfrag] = -1 # Read cannot extend beyond sink
            
            if j == len(ranges)-1: # Ending block
                if strand == 1 and e_tag: # Right position is a 3' end
                    tr = self.end_of_cluster(r, weight, self.end_ranges[1], self.end_extend, False)
                    if tr >= 0:
                        membership[sink_plus] = 1 # Add t+ to the membership table
                        r = tr
                        rfrag = self.frag_by_pos[r-1]
                        membership[(rfrag+1):(width-4)] = -1 # Read cannot extend beyond sink
                elif strand == -1 and s_tag: # Right position is a 5' end
                    tr = self.end_of_cluster(r, weight, self.end_ranges[2], self.end_extend*capped, capped)
                    if tr >= 0:
                        membership[source_minus] = 1 # Add s- to the membership table
                        r = tr
                        rfrag = self.frag_by_pos[r-1]
                        membership[(rfrag+1):(width-4)] = -1 # Read cannot extend beyond source
            
            if lfrag > rfrag: # Reassignment of ends caused lfrag and rfrag to be out of order
                if len(ranges) > 1:
                    if j == 0 and splice[j]: # The right border is a splice junction, structural violation
                        membership[:] = -1 # Read is a violation, remove
                        break
                    elif j == len(ranges)-1 and splice[j-1]: # The left border is a splice junction, structural violation
                        membership[:] = -1 # Read is a violation, remove
                        break
                
                if (s_tag and strand == 1) or (e_tag and strand == -1): # The left border was updated
                    rfrag = lfrag
                else: # The right border was updated
                    lfrag = rfrag
            
            membership[lfrag:(rfrag+1)] = 1 # Add all covered frags to the membership table
            if j > 0: # Processing a downstream block
                if splice[j-1]: # The gap between this block and the last was a splice junction
                    # Check that this junction is in the list of splice junctions
                    # If it was filtered out, this read is invalid and should be removed.
                    span = (self.frags[last_rfrag][1], self.frags[lfrag][0])
                    junction_hash = self.span_to_string(span)
                    if strand == 1 and junction_hash not in self.J_plus:
                        membership[:] = -1 # Read contains a filtered junction, remove
                        break
                    elif strand == -1 and junction_hash not in self.J_minus:
                        membership[:] = -1 # Read contains a filtered junction, remove
                        break
                    else:
                        membership[(last_rfrag+1):lfrag] = -1 # All frags in the intron are incompatible
                else: # The gap is an unspecified gap
                    intervening_junctions = self.junctions_between(self.frags[last_rfrag][1], self.frags[lfrag][0], strand)
                    if len(intervening_junctions) == 0: # Fill in the intervening gap if no junctions exist in the range
                        membership[(last_rfrag+1):lfrag] = 1
                    else: # Fill in as much information as possible based on filtered junctions and membership
                        spans = sorted([self.string_to_span(junction) for junction in intervening_junctions])
                        splice_sites = [site for span in spans for site in span] # 'unlist' the splice sites in the order they appear
                        membership[(last_rfrag+1):self.frag_by_pos[min(splice_sites)]] = 1 # Fill in up to the first junction boundary
                        membership[self.frag_by_pos[max(splice_sites)]:rfrag] = 1 # Fill in after the last junction boundary
                        plus_junctions = any([junction in self.J_plus.keys() for junction in intervening_junctions])
                        minus_junctions = any([junction in self.J_minus.keys() for junction in intervening_junctions])
                        if plus_junctions and minus_junctions:
                            break
                        elif plus_junctions:
                            jstrand = 1
                        else:
                            jstrand = -1
                        
                        sorted_splice_sites = sorted(set(splice_sites))
                        junctions_are_linear = splice_sites == sorted_splice_sites # Any nesting will cause this evaluation to be false
                        if junctions_are_linear: # No overlapping splice junctions
                            membership[(last_rfrag+1):lfrag] = 1
                            for span in spans:
                                skipped_frags = list(range(self.frag_by_pos[span[0]], self.frag_by_pos[span[1]]))
                                if not np.all(self.discard_frags[[strand>=0, strand<=0],:][:,skipped_frags]): # Two alternative paths exist through this intron (spliced and unspliced)
                                    for skipped in skipped_frags:
                                        membership[skipped] = 0
                                else: # Only the spliced path exists
                                    for skipped in skipped_frags:
                                        membership[skipped] = -1
                                    
                                    if strand == 0:
                                        strand = jstrand
                                        if strand == 1:
                                            membership[source_minus] = -1 # Read cannot have minus-stranded features
                                            membership[sink_minus] = -1 # Read cannot have minus-stranded features
                                        elif strand == -1:
                                            membership[source_plus] = -1 # Read cannot have plus-stranded features
                                            membership[sink_plus] = -1 # Read cannot have plus-stranded features
            
            last_rfrag = rfrag
        
        # Apply intron filtering to (a) restrict the read to one strand or (b) remove it entirely
        members = membership[:-4] == 1
        if not np.any(members):
            membership[:] = -1
        elif np.any(self.discard_frags[[strand>=0, strand<=0],:][:,members]):
            if strand != 0:
                membership[:] = -1 # Read contains a discarded frag, remove
            else: # The read may only be possible in one stranded orientation (or it may be impossible given discard_frags)
                if not np.any(self.discard_frags[0,members]): # Plus strand is still viable
                    membership[source_minus] = -1 # Read cannot have minus-stranded features
                    membership[sink_minus] = -1 # Read cannot have minus-stranded features
                elif not np.any(self.discard_frags[1,members]): # Minus strand is still viable
                    membership[source_plus] = -1 # Read cannot have minus-stranded features
                    membership[sink_plus] = -1 # Read cannot have minus-stranded features
                else: # Neither strand is viable, remove
                    membership[:] = -1

        return membership
    
    cpdef list junctions_between(self, int lpos, int rpos, char strand):
        """Returns a list of junction hashes that fall between the two specified positions
        on the specified strand (either strand if strand==0)."""
        cdef list candidates, junctions_in_range
        cdef str junction_hash
        cdef (int, int) span
        if strand == 1:
            candidates = list(self.J_plus.keys())
        elif strand == -1:
            candidates = list(self.J_minus.keys())
        else:
            candidates = list(self.J_plus.keys()) + list(self.J_minus.keys())
        
        junctions_in_range = list()
        for junction_hash in candidates:
            span = self.string_to_span(junction_hash)
            if span[0] >= lpos and span[1] <= rpos:
                junctions_in_range.append(junction_hash)
        
        return junctions_in_range
    
    cpdef void filter_members_by_strand(self):
        """If a read is in a region with >1-minimum_proportion coverage
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
    
    cpdef np.ndarray apply_intron_filter(self, float threshold=1):
        """Returns a 2-row array of bools marking non-viable frags for plus and minus paths
        if they are inferred to belong to an unprocessed transcript,
        i.e. a retained intron, run-on transcription downstream of a 3' end, or
        transcriptional noise upstream of a 5' end. The argument 'intron_filter' is used here."""
        cdef EndRange endrange
        cdef np.ndarray remove_plus, remove_minus, overlappers, flowthrough, terminal, sb, fills_intron, cov, junction_membership, discard_frags
        cdef list stranded_branches, nonterminal_frags
        cdef int strand, number_of_frags, frag, l, r, maxgap
        cdef str junction
        cdef (int, int) span
        # Check flowthrough across all starts/ends
        # Use self.frag_strand_ratios to assign non-stranded reads to strand-specific flowthroughs
        number_of_frags = len(self.frags)
        discard_frags = np.zeros((2,number_of_frags), dtype=np.bool)
        maxgap = self.extend
        nonterminal_frags = [(frag,span) for frag,span in enumerate(self.frags) if span[0] not in self.SPbp|self.EMbp and span[1] not in self.SMbp|self.EPbp]
        for frag,span in nonterminal_frags:
            if not passes_threshold(self.depth[span[0]:span[1]], maxgap, threshold):
                discard_frags[:,frag] = True
            # else:
            #     if not passes_threshold(self.cov_plus[self.frags[frag][0]:self.frags[frag][1]], maxgap, threshold):
            #         discard_frags[0,frag] = True
                
            #     if not passes_threshold(self.cov_minus[self.frags[frag][0]:self.frags[frag][1]], maxgap, threshold):
            #         discard_frags[1,frag] = True
        
        for endtype in range(4):
            if endtype < 2:
                cov = self.cov_plus
            else:
                cov = self.cov_minus
            
            for endrange in self.end_ranges[endtype]:
                frag = self.frag_by_pos[endrange.terminal - (endtype in [1,2])]
                if endtype == 0 and frag > 0: # S+, check plus-stranded left flowthrough
                    strand, overrun_frag = 1, frag-1
                    terminal_weight = cov[endrange.right-1]
                elif endtype == 1 and frag < number_of_frags-1: # E+, check plus-stranded right flowthrough
                    strand, overrun_frag = 1, frag+1
                    terminal_weight = cov[endrange.left]
                elif endtype == 2 and frag < number_of_frags-1: # S-, check minus-stranded right flowthrough
                    strand, overrun_frag = -1, frag+1
                    terminal_weight = cov[endrange.left]
                elif endtype == 3 and frag > 0: # E-, check minus-stranded left flowthrough
                    strand, overrun_frag = -1, frag-1
                    terminal_weight = cov[endrange.right-1]
                else:
                    continue
                
                overrun_cov = cov[self.frags[overrun_frag][0]:self.frags[overrun_frag][1]]
                flowthrough_weight = np.mean(overrun_cov)
                if np.any(self.depth[self.frags[overrun_frag][0]:self.frags[overrun_frag][1]] < threshold) or flowthrough_weight < self.intron_filter * terminal_weight:
                    if strand == 1:
                        discard_frags[0,overrun_frag] = True
                    else:
                        discard_frags[1,overrun_frag] = True
        
        # Check same-stranded intron retention
        strand = 1
        stranded_branches = list()
        for k in self.J_plus.keys():
            l,r = self.string_to_span(k)
            stranded_branches += [l,r]
        
        for i in [0,1]:
            for rng in self.end_ranges[i]:
                stranded_branches.append(rng.terminal)
        
        sb = np.unique(sorted(stranded_branches))
        for junction in self.J_plus.keys():
            l,r = self.string_to_span(junction)
            if not np.any(np.logical_and(sb>l, sb<r)): # Intron has no intervening stranded branchpoints
                lfrag = self.frag_by_pos[l]
                rfrag = self.frag_by_pos[r]
                cov_in_junction = self.depth[self.frags[lfrag][0]:self.frags[rfrag][0]]
                spanning_cov = np.max(np.append(self.depth[self.frags[lfrag][0]], self.depth[self.frags[rfrag][0]]))
                if np.any(cov_in_junction < threshold) or np.mean(cov_in_junction) < spanning_cov*self.intron_filter:
                    discard_frags[0, lfrag:rfrag] = True
        
        # Repeat the procedure for the minus strand
        strand = -1
        stranded_branches = list()
        for k in self.J_minus.keys():
            l,r = self.string_to_span(k)
            stranded_branches += [l,r]
        
        for i in [2,3]:
            for rng in self.end_ranges[i]:
                stranded_branches.append(rng.terminal)
        
        sb = np.unique(sorted(stranded_branches))
        for junction in self.J_minus.keys():
            l,r = self.string_to_span(junction)
            if not np.any(np.logical_and(sb>l, sb<r)): # Intron has no intervening stranded branchpoints
                lfrag = self.frag_by_pos[l]
                rfrag = self.frag_by_pos[r]
                cov_in_junction = self.depth[self.frags[lfrag][0]:self.frags[rfrag][0]]
                spanning_cov = np.max(np.append(self.depth[self.frags[lfrag][0]], self.depth[self.frags[rfrag][0]]))
                if np.any(cov_in_junction < threshold) or np.mean(cov_in_junction) < spanning_cov*self.intron_filter:
                    discard_frags[1, lfrag:rfrag] = True
        
        return discard_frags        
        
    cpdef np.ndarray get_competitors(self, int index):
        """Given an element index, return a list of all elements that 
        (1) share at least one member and 
        (2) are incompatible."""
        cdef np.ndarray incompatible, members, competitors
        incompatible =  np.where(self.overlap[:,index]==-1)[0]
        if len(incompatible) == 0:
            return incompatible
        
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
            new_reps = np.zeros(shape=reduced_membership.shape[0], dtype=np.float32)
            new_lengths = np.zeros(shape=reduced_membership.shape[0], dtype=np.int32)
            if not np.any(self.member_weights):
                member_weights_exists = False
            else:
                member_weights_exists = True
                new_member_weights = np.zeros(shape=(reduced_membership.shape[0],reduced_membership.shape[1]), dtype=np.float32)
            
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
            self.rep_array = new_reps[sorted_indices]
            if not member_weights_exists:
                self.member_weights = np.full((self.membership.shape[0],self.membership.shape[1]), np.sum(self.weight_array,axis=1,keepdims=True))
                self.member_weights[self.membership==0] = 0
                for i in range(self.membership.shape[0]):
                    self.member_weights[i,self.membership[i,:]==-1] = self.rep_array[i]
            else:
                self.member_weights = new_member_weights[sorted_indices,:]
            
            self.member_lengths = new_lengths[sorted_indices]
            self.strand_array = new_strands[sorted_indices]
        
        if not np.any(self.member_weights): # member_weights still uninitialized
            self.member_weights = np.full((self.membership.shape[0],self.membership.shape[1]), np.sum(self.weight_array,axis=1,keepdims=True))
            self.member_weights[self.membership==0] = 0
            for i in range(self.membership.shape[0])    :
                self.member_weights[i,self.membership[i,:]==-1] = self.rep_array[i]
    
    cpdef void filter_by_reps(self, float threshold=1):
        """Enforce that elements in the membership"""
        cdef np.ndarray keep
        if threshold > 1:
            keep = np.where(self.rep_array > threshold)[0]
            self.subset_elements(keep)
    
    cpdef void build_overlap_matrix(self):
        """Builds reduced-representation matrix where every read is (-1) incompatible,
        (1) overlapping, or (0) non-overlapping with every other read. Any read completely
        contained in one or more reads with more information are removed, with their
        weight distributed proportionally to the containers.
        """
        if self.ignore_ends:
            endless_matrix = remove_ends(self.membership)
            endless_info = get_information_content(endless_matrix)
            self.overlap = calculate_overlap_matrix(endless_matrix, endless_info, self.strand_array)
        else:
            self.overlap = calculate_overlap_matrix(self.membership[:,[-4,-1]+list(range(self.membership.shape[1]-4))+[-3,-2]], self.information_content, self.strand_array)
    
    cpdef void subset_elements(self, np.ndarray keep):
        if not self.membership is None: self.membership = self.membership[keep,:]
        if not self.overlap is None: self.overlap = self.overlap[keep,:][:,keep]
        if not self.weight_array is None: self.weight_array = self.weight_array[keep,:]
        if not self.member_weights is None: self.member_weights = self.member_weights[keep,:]
        if not self.rep_array is None: self.rep_array = self.rep_array[keep]
        if not self.strand_array is None: self.strand_array = self.strand_array[keep]
        if not self.information_content is None: self.information_content = self.information_content[keep]
        if not self.member_content is None: self.member_content = self.member_content[keep]
        if not self.member_lengths is None: self.member_lengths = self.member_lengths[keep]
        self.bases = np.sum(np.sum(self.weight_array, axis=1)*self.member_lengths)
        self.number_of_elements = len(keep)
    
    cpdef void collapse_chains(self):
        """Split the Overlap Matrix into a list of connected components.
        Each component can be assembled independently."""
        # cdef strandedComponents cc
        cdef simplifyDFS dfs
        cdef list simplified_indices
        cdef dict chains
        cdef np.ndarray indices
        dfs = simplifyDFS(self.overlap, np.argsort(self.information_content))
        simplified_indices = []
        chains = {}
        for i in range(dfs.vertices):
            chain = dfs.component[i]
            if chain in chains.keys():
                parent = chains[chain]
                self.merge_reads(i, parent)
            else:
                simplified_indices.append(i)
                chains[chain] = i
        
        indices = np.array(simplified_indices, dtype=np.int32)
        self.subset_elements(indices)
        self.overlap = calculate_overlap_matrix(self.membership[:,[-4,-1]+list(range(self.membership.shape[1]-4))+[-3,-2]], self.information_content, self.strand_array)
        # self.adj = {i:[] for i in range(self.overlap.shape[0])}
        # edge_locations = np.where(self.overlap >= 1)
        # for a,b in zip(edge_locations[0],edge_locations[1]):
        #     if a != b:
        #         self.adj[a].append(b)
        
        # self.exc = {i:set(np.where(self.overlap[i,:]==-1)[0]) for i in range(self.overlap.shape[0])}
        # cc = strandedComponents(self.adj, self.strand_array)
        # component_bool = np.zeros((self.membership.shape[0], cc.c), dtype=np.bool)
        # for c in cc.pc:
        #     indices = cc.component_plus==c
        #     prior_overlap = np.any(component_bool[indices,:c],axis=0)
        #     if not self.ignore_ends:
        #         if not np.all(np.sum(self.membership[indices,-4:-2]==1,axis=0)>0): # At least one start and end exists in the component
        #             continue
            
        #     component_bool[:,c] = indices
        #     if np.any(prior_overlap):
        #         prior_components = np.where(prior_overlap)[0]
        #         component_bool[:,c] += np.sum(component_bool[:,prior_components],axis=1)>0
        #         component_bool[:,prior_components] = False
        
        # for c in cc.mc:
        #     indices = cc.component_minus==c
        #     prior_overlap = np.any(component_bool[indices,:c],axis=0)
        #     if not self.ignore_ends:
        #         if not np.all(np.sum(self.membership[indices,-2:]==1,axis=0)>0): # At least one start and end exists in the component
        #             continue
            
        #     component_bool[:,c] = indices
        #     if np.any(prior_overlap):
        #         prior_components = np.where(prior_overlap)[0]
        #         component_bool[:,c] += np.sum(component_bool[:,prior_components],axis=1)>0
        #         component_bool[:,prior_components] = False
        
        # subproblems = []
        # for c in range(component_bool.shape[1]):
        #     if np.any(component_bool[:,c]):
        #         indices = np.where(component_bool[:,c])[0]
        #         dfs = simplifyDFS(self.overlap[indices,:][:,indices], np.argsort(self.information_content[indices]))
        #         simplified_indices = []
        #         chains = {}
        #         for i in range(dfs.vertices):
        #             chain = dfs.component[i]
        #             if chain in chains.keys():
        #                 parent = chains[chain]
        #                 self.merge_reads(indices[i], parent)
        #             else:
        #                 simplified_indices.append(indices[i])
        #                 chains[chain] = indices[i]
                
        #         indices = np.array(simplified_indices, dtype=np.int32)
        #         subproblems += [indices]
        
        # self.overlap = calculate_overlap_matrix(self.membership, self.information_content, self.strand_array)
        # return subproblems
    
    cpdef void build_graph(self, reduce=True):
        """Constructs one or more graphs from connection values (ones) in the overlap matrix.
        Each graph is an _element_graph.ElementGraph() object with a built-in assembly method.
        """
        if reduce: # Split graph into connected components and solve each on its own
            self.collapse_chains()
        
        self.graph = ElementGraph(self.overlap, self.membership, self.weight_array, self.member_weights, self.strand_array, self.frag_len, self.naive, dead_end_penalty=self.dead_end_penalty, ignore_ends=self.ignore_ends, intron_filter=self.intron_filter)
    
    cpdef void assemble_transcripts(self):
        if self.graph is not None:
            self.graph.assemble(self.minimum_proportion)
            counter = 1
            for path in self.graph.paths:
                self.transcripts += self.convert_path(path, counter)
                counter += 1
            
            self.add_transcript_attributes()
    
    cpdef void trim_transcript_ends(self, transcript):
        """(for endless transcripts only) Cuts 5' and 3' down
        to the position with the largest relative delta."""
        cdef np.ndarray d, dd, posdelta, negdelta, dist_from_end
        cdef bint update_left, update_right
        cdef int l, r, maxdelta
        update_left, update_right = False, False
        if transcript.strand == 0:
            update_left = True
            update_right = True
        
        if (transcript.strand == 1 and not transcript.s_tag) or (transcript.strand == -1 and not transcript.e_tag):
            update_left = True
        
        if (transcript.strand == 1 and not transcript.e_tag) or (transcript.strand == -1 and not transcript.s_tag):
            update_right = True
        
        if update_left:
            l = transcript.ranges[0][0] - self.leftmost
            r = transcript.ranges[0][1] - self.leftmost
            d = self.depth[l:r]
            dd = np.diff(d)
            posdeltas = np.where(np.logical_and(dd > 0,d[1:]<np.mean(d)))[0]
            if len(posdeltas) > 0:
                dist_from_end = 1 - posdeltas / (r-l)
                maxdelta = posdeltas[np.argsort(dist_from_end * -np.power(dd[posdeltas],2)/d[posdeltas])[0]]
                transcript.ranges[0] = (transcript.ranges[0][0]+maxdelta, transcript.ranges[0][1])
        
        if update_right:
            l = transcript.ranges[-1][0] - self.leftmost
            r = transcript.ranges[-1][1] - self.leftmost
            d = self.depth[l:r]
            dd = np.diff(d)
            negdeltas = np.where(np.logical_and(dd < 0, d[:-1] < np.mean(d)))[0]
            if len(negdeltas) > 0:
                dist_from_end = negdeltas / (r-l)
                maxdelta = negdeltas[np.argsort(dist_from_end * -np.power(dd[negdeltas],2)/d[negdeltas-1])[0]]
                transcript.ranges[-1] = (transcript.ranges[-1][0], transcript.ranges[-1][0]+maxdelta)
        
        return

    cpdef void add_transcript_attributes(self):
        """Populate the new read objects with diagnostic information
        to store in the GTF attributes column."""
        cdef:
            int first, last, s_pos, e_pos
            dict S_info, E_info
            list S_ranges, E_ranges
            EndRange S, E
            (int, int) span
        
        for T in self.transcripts:
            T.attributes['length'] = T.get_length()
            S_info = {'S.reads':0, 'S.capped':0, 'S.left':0, 'S.right':0}
            E_info = {'E.reads':0, 'E.left':0, 'E.right':0}
            
            if T.strand != 0:
                if T.strand == 1:
                    s_pos = first = T.span[0] - self.leftmost
                    e_pos = last = T.span[1] - self.leftmost
                    S_ranges = self.end_ranges[0]
                    E_ranges = self.end_ranges[1]
                else:
                    s_pos = first = T.span[1] - self.leftmost
                    e_pos = last = T.span[0] - self.leftmost
                    S_ranges = self.end_ranges[2]
                    E_ranges = self.end_ranges[3]
                
                if T.s_tag:
                    S = self.get_end_cluster(first, 0, S_ranges, self.end_extend)
                    if S is not self.nullRange:
                        s_pos = S.peak
                        span = S.span()
                        S_info['S.reads'] = round(sum(S.positions.values()),1)
                        S_info['S.capped'] = round(S.capped,1)
                        S_info['S.left'] = span[0] + self.leftmost
                        S_info['S.right'] = span[1] + self.leftmost + 1
                        if s_pos != first: # S pos was replaced
                            if T.strand == 1:
                                T.ranges[0] = (s_pos + self.leftmost, T.ranges[0][1])
                            else:
                                T.ranges[-1] = (T.ranges[-1][0], s_pos + self.leftmost + 1)
                
                if T.e_tag:
                    E = self.get_end_cluster(last, 0, E_ranges, self.end_extend)
                    if E is not self.nullRange:
                        e_pos = E.peak
                        span = E.span()
                        E_info['E.reads'] = round(sum(E.positions.values()),1)
                        E_info['E.left'] = span[0] + self.leftmost
                        E_info['E.right'] = span[1] + self.leftmost + 1
                        if e_pos != last: # S pos was replaced
                            if T.strand == 1:
                                T.ranges[-1] = (T.ranges[-1][0], e_pos + self.leftmost + 1)
                            else:
                                T.ranges[0] = (e_pos + self.leftmost, T.ranges[0][1])
            
            T.attributes.update(S_info)
            T.attributes.update(E_info)
            if T.attributes['S.capped'] > 0 and T.attributes['S.capped'] >= T.attributes['S.reads']*.1:
                T.capped = True
    
    cpdef void merge_reads(self, int child_index, int parent_index):
        """Combines the information of two read elements in the object."""
        cdef char p_out, c_out, p_in, c_in, s
        cdef (char, char) o
        cdef Py_ssize_t i
        cdef int combined_length
        # Parent gains all information of child
        for i in range(self.membership.shape[1]): # Iterate over columns of the membership table
            p = self.membership[parent_index,i]
            c = self.membership[child_index,i]  
            if p == 0: # If no info in parent, overwrite with child
                self.membership[parent_index, i] = c
            elif c != 0 and p != c: # Conflicting membership information
                raise Exception('Incompatible read pair: {}, {}'.format(child_index, parent_index))
        
        self.information_content[parent_index] = np.sum(np.abs(self.membership[parent_index,:]))
        s = self.strand_array[parent_index]
        if s == 0:
            self.strand_array[parent_index] = self.strand_array[child_index]
        
        combined_length = np.sum(self.frag_len[self.membership[parent_index,:]==1])
        self.weight_array[parent_index,:] = (self.weight_array[child_index,:]*self.member_lengths[child_index] + self.weight_array[parent_index,:]*self.member_lengths[parent_index])/combined_length
        self.member_weights[parent_index,:] += self.member_weights[child_index,:]
        self.member_lengths[parent_index] = combined_length
        self.weight_array[child_index,:] = 0
        self.member_weights[child_index,:] = 0
    
    cpdef list convert_path(self, element, transcript_number):
        """Prints a representation of an ElementGraph Element object
        by converting it first to an RNAseqMapping object."""
        cdef int chrom, source, N, m, n, l, r, last_member
        cdef char strand
        cdef bint s_tag, e_tag, capped, gap_is_splice
        cdef list ranges, splice, members, output
        cdef (int, int) frag, exon
        cdef set nonmembers
        cdef str gene_id, transcript_id, junction_hash
        cdef dict junctions
        cdef EndRange S, E
        output = []
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
                    S = self.get_end_cluster(l, 0, self.end_ranges[0], self.end_extend)
                    if S is not self.nullRange:
                        l = S.peak
                elif element.e_tag and element.strand == -1:
                    E = self.get_end_cluster(l, 0, self.end_ranges[3], self.end_extend)
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
            S = self.get_end_cluster(r-1, 0, self.end_ranges[2], self.end_extend)
            if S is not self.nullRange:
                r = S.peak+1
        elif element.e_tag and element.strand == 1:
            E = self.get_end_cluster(r-1, 0, self.end_ranges[1], self.end_extend)
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
        readObject.attributes['cov'] = round(element.bases/element.length, 2)
        readObject.attributes['bases'] = round(element.bases, 2)
        readObject.coverage = readObject.attributes['cov']
        if self.allow_incomplete and not readObject.complete:
            self.trim_transcript_ends(readObject)
        
        return [readObject]


cdef class simplifyDFS():
    cdef public dict O, X, CO, CX
    cdef public int vertices, c
    cdef public np.ndarray visited, component, pre, post
    def __init__(self, overlap_matrix, search_order):
        self.O = self.getOutgroups(overlap_matrix)
        self.X = {i:set(np.where(overlap_matrix[i,:]==-1)[0]) for i in range(overlap_matrix.shape[0])}
        self.CO = {} # Outgroups of each component
        self.CX = {} # Exclusions of each component
        self.vertices = len(self.O.keys())
        self.visited = np.zeros(self.vertices, dtype=np.bool)
        self.component = np.full(self.vertices, -1, dtype=np.int32)
        self.pre = np.zeros(self.vertices, dtype=np.int32)
        self.post = np.zeros(self.vertices, dtype=np.int32)
        self.c = 0
        clock = 0
        for v in search_order:
            if not self.visited[v]:
                clock = self.Explore(v, clock)
    
    cdef dict getOutgroups(self, np.ndarray overlap_matrix):
        cdef tuple edge_locations
        cdef dict O
        O = {i:[] for i in range(overlap_matrix.shape[0])} # Outgroups of each element
        edge_locations = np.where(overlap_matrix >= 1)
        for a,b in zip(edge_locations[0],edge_locations[1]):
            if a != b:
                O[a].append(b)
        
        return O
    
    cdef int Previsit(self, int v, int clock):
        self.pre[v] = clock
        return clock + 1
    
    cdef void makeComponent(self, int v):
        self.c += 1
        self.component[v] = self.c
        self.CO[self.c] = set(self.O[v]+[v])
        self.CX[self.c] = self.X[v]
    
    cdef int Postvisit(self, int v, int clock):
        self.post[v] = clock
        outgroups = np.unique(self.component[self.O[v]])
        outgroups = outgroups[outgroups!=-1]
        for outgroup in outgroups: # Check if v can be added to the component of any of its outgroups
            if set(self.O[v]).issubset(self.CO[outgroup]) and self.X[v] == self.CX[outgroup]:
                self.component[v] = outgroup
                self.CO[outgroup].add(v)
                return clock + 1
        
        # If v is compatible with no outgroups, add a new component
        self.makeComponent(v)
        return clock + 1
    
    cdef int Explore(self, int v, int clock):
        cdef int w
        self.visited[v] = True
        clock = self.Previsit(v, clock)
        for w in self.O[v]:
            # print('{}->{} ({})'.format(v,w,self.visited[w]))
            if not self.visited[w]:
                # print("Visiting {}".format(w))
                clock = self.Explore(w, clock)
        
        # print("Closing {}".format(v))
        clock = self.Postvisit(v, clock)
        # print(self.component)
        return clock

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


cpdef (char,char) get_overlap(np.ndarray[char, ndim=1] members_a, np.ndarray[char, ndim=1] members_b, int info_a, int info_b):
    """Returns the a->b and b->a overlap relationship between two reads"""
    cdef int ia, ib, shared, a_to_b, b_to_a
    cdef (bint, bint, bint, bint) info_buffer
    cdef bint overlapping
    cdef char horiz, vert
    info_buffer = (False, False, False, False)
    shared, a_to_b, b_to_a = 0,0,0
    overlapping, in_a, in_b = False, False, False
    for i in range(len(members_a)):
        ia = members_a[i]
        ib = members_b[i]
        if (ia == 1 and ib == -1) or (ia == -1 and ib == 1): # Incompatibility found
            return (-1, -1)
        
        # If not incompatible, then they either share or do not share membership
        shared += ia == ib and ia != 0 # Shared information (inclusion or exclusion)
        overlapping = overlapping or ia + ib == 2 # At least one member is shared
        info_buffer = (ia!=0, info_buffer[0], ib!=0, info_buffer[2])
        a_to_b += in_a and info_buffer == (False, True, True, True)
        b_to_a += in_b and info_buffer == (True, True, False, True)
        in_a = ia!=0 and (in_a or ia==1)
        in_b = ib!=0 and (in_b or ib==1)
    
    if shared <= 0:
        return (0, 0)
    
    if shared == info_a:
        horiz = 2
    elif shared == info_b:
        vert = 2
    else:
        horiz = int(overlapping and a_to_b > 0)
        vert = int(overlapping and b_to_a > 0)

    return (horiz, vert) 

cpdef np.ndarray calculate_overlap_matrix(np.ndarray[char, ndim=2] membership_matrix, np.ndarray information_content, np.ndarray strand_array):
    """Given a matrix of membership values (1, 0, or -1; see self.build_membership_matrix),
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
    cdef char sa, sb, horiz, vert
    cdef Py_ssize_t a, b, maxlen, number_of_reads, number_of_frags
    
    number_of_reads = membership_matrix.shape[0]
    number_of_frags = membership_matrix.shape[1]
    if len(information_content) != number_of_reads:
        information_content = get_information_content(membership_matrix)
    
    cdef np.ndarray[char, ndim=2] overlap_matrix = np.zeros((number_of_reads, number_of_reads), dtype=np.int8) # Container for overlap information
    cdef char [:] STRAND_ARRAY = strand_array
    cdef int [:] INFO = information_content
    cdef char [:, :] COMPATIBILITY = overlap_matrix
    maxlen = COMPATIBILITY.shape[0]
    for a in range(maxlen):
        for b in range(a,maxlen):
            if b == a: # A read necessarily contains itself
                COMPATIBILITY[a,b] = 2
                continue
            
            sa = STRAND_ARRAY[a]
            sb = STRAND_ARRAY[b]
            if (sa == 1 and sb == -1) or (sa == -1 and sb == 1):
                COMPATIBILITY[a,b] = -1
                COMPATIBILITY[b,a] = -1
                continue
            
            horiz, vert = get_overlap(membership_matrix[a], membership_matrix[b], INFO[a], INFO[b])
            COMPATIBILITY[a,b] = horiz
            COMPATIBILITY[b,a] = vert
    
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

def sum_subset(mask, array_to_mask):
    """Given an array and a boolean mask, return
    the sum of array values at which mask was True"""
    return array_to_mask[mask].sum()

