#cython: language_level=3
import numpy as np
cimport numpy as np
import re
import copy
import json
from bookend.core.cython_utils._element_graph import ElementGraph
from bookend.core.cython_utils._rnaseq_utils import RNAseqMapping, ELdata, range_of_reads, build_depth_matrix
from collections import deque, Counter
import cython
import time

bp_typeorder = {'N':-1, 'E':0, 'A':1, 'D':2, 'S':3} # Sort order for branchpoint types

cdef class Locus:
    cdef public int chrom, leftmost, rightmost, extend, end_extend, number_of_elements, min_overhang, chunk_number
    cdef public bint naive, infer_starts, infer_ends, use_attributes
    cdef public tuple reads, frags
    cdef public float weight, minimum_proportion, cap_bonus, novelty_ratio, mean_read_length, intron_filter
    cdef public dict adj, bp_lookup, J_plus, J_minus
    cdef public list transcripts, traceback
    cdef public object BP, graph
    cdef public np.ndarray depth_matrix, strandscaled, cov_plus, cov_minus, read_lengths, frag_len, frag_by_pos, strand_array, weight_array, rep_array, membership, overlap, information_content, member_content
    def __init__(self, chrom, chunk_number, list_of_reads, extend=0, end_extend=100, min_overhang=3, reduce=True, minimum_proportion=0.02, cap_bonus=5, novelty_ratio=1, complete=False, verbose=False, naive=True, intron_filter=0.15, infer_starts=False, infer_ends=False, use_attributes=False):
        self.transcripts = []
        self.traceback = []
        self.bp_lookup = {}
        self.chunk_number = chunk_number
        self.naive = naive
        self.minimum_proportion = minimum_proportion
        self.intron_filter = intron_filter
        self.min_overhang = min_overhang
        self.chrom = chrom
        self.novelty_ratio = novelty_ratio
        self.cap_bonus = cap_bonus
        self.infer_starts = infer_starts
        self.infer_ends = infer_ends
        self.use_attributes = use_attributes
        if len(list_of_reads) > 0:
            self.leftmost, self.rightmost = range_of_reads(list_of_reads)
            self.reads = tuple(list_of_reads) # Cannot be mutated
            self.read_lengths = np.array([r.get_length() for r in self.reads])
            self.mean_read_length = np.mean(self.read_lengths)
            self.weight = float(0)
            self.extend = extend
            self.end_extend = end_extend
            self.depth_matrix, self.J_plus, self.J_minus = build_depth_matrix(self.leftmost, self.rightmost, list_of_reads, self.cap_bonus, self.use_attributes)
            self.generate_branchpoints()
            # if type(self) is AnnotationLocus:
            #     self.traceback = [set([i]) for i in range(len(self.reads))]
            #     self.build_membership_matrix(0)
            #     self.filter_by_reps(self.minreps)
            #     if self.membership.shape[0] > 0:
            #         self.build_overlap_matrix(reduce=False, ignore_ends=True)
            # else:
            #     self.build_membership_matrix()
            #     self.denoise(verbose)
            #     self.filter_run_on_frags()
            #     if self.membership.shape[0] > 0:
            #         self.build_overlap_matrix(reduce)
            #         self.build_graph(reduce)
    
    def __len__(self):
        return self.rightmost - self.leftmost
    
    def __repr__(self):
        symbols = {-1:'-', 0:' ', 1:'+', 2:'^'}
        summary_string = '<{} ({})>\n'.format(str(type(self)).split("'")[-2], self.number_of_elements)
        for l in range(self.number_of_elements):
            members = ''.join([symbols[i] for i in self.membership[l,:]])
            overlap = ''.join([symbols[i] for i in self.overlap[l,:]])
            summary_string += '  |{}|\t|{}|\n'.format(members, overlap)
        
        return summary_string

    cpdef generate_branchpoints(self):
        """Estimate strand-specific coverage, then use it to generate an ordered array of
        positions where branching could occur in the overlap graph."""
        cdef:
            np.ndarray strandratio, covstranded, strandedpositions, pos, vals, value_order
            Py_ssize_t Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn, covrow, i
            float threshold
        
        Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn = range(11)
        strandratio = np.full(self.depth_matrix.shape[1], 0.5, dtype=np.float32)
        covstranded = np.sum(self.depth_matrix[(covp,covm),:],axis=0)
        strandedpositions = np.where(covstranded > 0)[0]
        strandratio[strandedpositions] = self.depth_matrix[covp,strandedpositions]/covstranded[strandedpositions]
        self.cov_plus = self.depth_matrix[covp,] + self.depth_matrix[covn,]*strandratio
        self.cov_minus = self.depth_matrix[covm,] + self.depth_matrix[covn,]*(1-strandratio)
        end_ranges = dict()
        for endtype in [Sp, Ep, Sm, Em]:
            pos = np.where(self.depth_matrix[endtype,]>0)[0]
            if endtype in [Sp, Ep]:
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_plus[pos]
            else:
                vals = np.power(self.depth_matrix[endtype, pos],2)/self.cov_minus[pos]
            
            end_ranges[endtype] = self.make_end_ranges(pos, vals)
        
        self.branchpoints = end_ranges

    cpdef list make_end_ranges(self, np.ndarray pos, np.ndarray vals):
        """Returns a list of tuples that (1) filters low-signal positions
        and (2) clusters high-signal positions within self.end_extend.
        Returns list of (l,r) tuples demarking the edges of clustered end positions."""
        cdef:
            np.ndarray value_order, filtered_pos
            int p, ia
            float cumulative, threshold
            list ranges
        
        value_order = np.argsort(-vals)
        threshold = (1 - self.minimum_proportion) * np.sum(vals)
        cumulative = 0
        i = 0
        while cumulative < threshold:
            cumulative += vals[value_order[i]]
            i += 1
        
        filtered_pos = np.sort(pos[value_order[:i]])
        if filtered_pos.shape[0] == 0:
            return []
        elif filtered_pos.shape[0] == 1:
            p = filtered_pos[0]
            return [(p, p)]
        else:
            p = filtered_pos[0]
            ranges = [(p, p)]
            for p in filtered_pos[1:]:
                if p - self.end_extend <= ranges[-1][1]:
                    ranges[-1] = (ranges[-1][0], p)
                else:
                    ranges += [(p, p)]
        
        return ranges
    
    cdef str span_to_string(self, (int, int) span):
        """Converts a tuple of two ints to a string connected by ':'"""
        return '{}:{}'.format(span[0], span[1])
    
    cdef str string_to_span(self, str string):
        """Converts a string from span_to_string() back into a span"""
        cdef list splitstring = string.split(':')
        return (int(splitstring[0]), int(splitstring[1]))
    
    cpdef build_membership_matrix(self, threshold=1):
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
        cdef int last_rfrag, l, r, l_overhang, r_overhang, lfrag, rfrag, pos, locus_length
        cdef np.ndarray membership, strand_array, weight_array
        cdef char s
        cdef list temp_frags
        cdef (int, int) block, span
        cdef str junction_hash
        cdef set positions = set()
        
        bp_positions = sorted(list(positions))
        temp_frags = []
        for i in range(len(bp_positions)-1):
            block = (bp_positions[i], bp_positions[i+1])
            temp_frags.append(block)
        
        self.frags = tuple(temp_frags) # Immutable after init, makes a vertex between all pairs of nonterminal branchpoint positions
        self.frag_len = np.array([b-a for a,b in self.frags]+[int(self.mean_read_length*.5)]*4, dtype=np.int32)
        # self.frag_len = np.array([b-a for a,b in self.frags]+[0,0,0,0], dtype=np.int32)
        self.frag_by_pos = np.full(shape=len(self), fill_value=-1, dtype=np.int32)
        
        for i in range(len(self.frags)):
            if i == 0:
                a = i
            else:
                a = self.frags[i][0] - self.leftmost
            
            if i == len(self.frags) - 1:
                b = len(self.frag_by_pos)
            else:
                b = self.frags[i][1] - self.leftmost
                
            self.frag_by_pos[a:b] = i
        
        number_of_reads = len(self.reads)
        number_of_frags = len(self.frags)
        self.rep_array = np.ones(number_of_reads, dtype=np.int32)
        membership = np.zeros((number_of_reads, number_of_frags+4), dtype=np.int8) # Container for membership of each vertex in each read (-1 False, 1 True, 0 Unmeasured)
        strand_array = np.zeros(number_of_reads, dtype=np.int8) # Container for strandedness of each read (-1 minus, 1 plus, 0 nonstranded)
        weight_array = np.array([read.weight for read in self.reads], dtype=np.float32)
        cdef char [:, :] MEMBERSHIP = membership
        source_plus, sink_plus, source_minus, sink_minus = range(number_of_frags, number_of_frags+4)
        locus_length = len(self.frag_by_pos)
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
                    if self.frag_len[lfrag] > self.min_overhang:
                        if l+self.min_overhang < locus_length:
                            l_overhang = self.frag_by_pos[l+self.min_overhang-1]
                            if l_overhang != lfrag: # The left overhang is too short
                                lfrag = l_overhang
                    
                    if self.frag_len[rfrag] > self.min_overhang:
                        if r-self.min_overhang >= 0:
                            r_overhang = self.frag_by_pos[r-self.min_overhang]
                            if r_overhang != rfrag: # The left overhang is too short
                                rfrag = r_overhang

                if j == 0: # Starting block
                    if s == 1 and read.s_tag: # Left position is a 5' end
                        if l+self.leftmost in self.BP.S_plus.keys():
                            MEMBERSHIP[i, source_plus] = 1 # Add s+ to the membership table
                            l = self.BP.branchpoints[self.BP.S_plus[l+self.leftmost]].pos - self.leftmost
                            lfrag = self.frag_by_pos[l]
                            MEMBERSHIP[i, 0:lfrag] = -1 # Read cannot extend beyond source
                    elif s == -1 and read.e_tag: # Left position is a 3' end
                        if l+self.leftmost in self.BP.E_minus.keys():
                            MEMBERSHIP[i, sink_minus] = 1 # Add t- to the membership table
                            l = self.BP.branchpoints[self.BP.E_minus[l+self.leftmost]].pos - self.leftmost
                            lfrag = self.frag_by_pos[l]
                            MEMBERSHIP[i, 0:lfrag] = -1 # Read cannot extend beyond sink

                if j == len(read.ranges)-1: # Ending block
                    if s == 1 and read.e_tag: # Right position is a 3' end
                        if r+self.leftmost in self.BP.E_plus.keys():
                            MEMBERSHIP[i, sink_plus] = 1 # Add t+ to the membership table
                            r = self.BP.branchpoints[self.BP.E_plus[r+self.leftmost]].pos - self.leftmost
                            rfrag = self.frag_by_pos[r-1]
                            MEMBERSHIP[i, (rfrag+1):number_of_frags] = -1 # Read cannot extend beyond sink
                    elif s == -1 and read.s_tag: # Right position is a 5' end
                        if r+self.leftmost in self.BP.S_minus.keys():
                            MEMBERSHIP[i, source_minus] = 1 # Add s- to the membership table
                            r = self.BP.branchpoints[self.BP.S_minus[r+self.leftmost]].pos - self.leftmost
                            rfrag = self.frag_by_pos[r-1]
                            MEMBERSHIP[i, (rfrag+1):number_of_frags] = -1 # Read cannot extend beyond source

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
                        if s == 1 and junction_hash not in self.BP.J_plus:
                            MEMBERSHIP[i,:] = -1 # Read contains a filtered junction, remove
                            break
                        elif s == -1 and junction_hash not in self.BP.J_minus:
                            MEMBERSHIP[i,:] = -1 # Read contains a filtered junction, remove
                            break
                        
                        MEMBERSHIP[i, (last_rfrag+1):lfrag] = -1 # All frags in the intron are incompatible

                last_rfrag = rfrag

        if threshold > 0:
            for i in range(len(self.frags)):
                l,r = self.frags[i]
                frag_depth = self.BP.depth[l-self.leftmost:r-self.leftmost]
                if not passes_threshold(frag_depth, self.extend, threshold): # This frag has too large of a gap
                    MEMBERSHIP[:,i] = -1
        
        self.membership = membership
        self.weight_array = weight_array
        self.strand_array = strand_array
        self.reduce_membership()
        self.weight = np.sum(self.weight_array)
        self.number_of_elements = self.membership.shape[0]
        self.information_content = get_information_content(self.membership)
        self.member_content = get_member_content(self.membership)

    cpdef void filter_run_on_frags(self):
        """Iterates over frags looking for those putatively connecting the ends
        of two transcripts together. Applies the stringent filter (intron_filter)
        to frags with the following pairs:
            E>E<, S<S> - converging, diverging
            E>S>, S<E< - consecutive
        Coverage of this frag must exceed intron_filter*coverage of the heavier side.
        If not, membership of this frag is set to all -1's.
        Reads that lose all their members are thrown out.
        """
        cdef:
            Py_ssize_t frag_index, bp_index, number_of_frags, number_of_bps, i
            int l, r
            tuple valid_end_pairs
            str end_pair
            bint removed_a_frag
            np.ndarray has_frag
        
        removed_a_frag = False
        #TODO: Check upstream of EVERY S and downstream of EVERY E
        valid_end_pairs = ('S<S>', 'E>E<', 'E>S>', 'S<E<')
        frag_index = bp_index = 0
        number_of_frags = len(self.frags)
        number_of_bps = len(self.BP.branchpoints)
        for frag_index in range(1,number_of_frags-1):
            l, r = self.frags[frag_index]
            lbp = self.bp_lookup[l]
            rbp = self.bp_lookup[r]
            end_pair = lbp.branchtype
            end_pair += '>' if lbp.strand == 1 else '<'
            end_pair += rbp.branchtype
            end_pair += '>' if rbp.strand == 1 else '<'
            if end_pair in valid_end_pairs:
                ll,lr = self.frags[frag_index-1]
                rl,rr = self.frags[frag_index+1]
                weight_with_frag = np.mean(self.BP.depth[l-self.leftmost:r-self.leftmost])
                flank_left = np.mean(self.BP.depth[ll-self.leftmost:lr-self.leftmost])
                flank_right = np.mean(self.BP.depth[rl-self.leftmost:rr-self.leftmost])
                if weight_with_frag < max(flank_left, flank_right)*self.intron_filter:
                    has_frag = np.where(self.membership[:,frag_index]==1)[0]
                    self.membership[has_frag,:] = -1
                    self.membership[:,frag_index] = -1
                    removed_a_frag = True
        
        if removed_a_frag:
            keep = np.where(np.sum(self.membership==1,axis=1)>0)[0]
            self.membership = self.membership[keep,:]
            self.weight_array = self.weight_array[keep]
            self.strand_array = self.strand_array[keep]
            self.member_content = self.member_content[keep]
            self.weight = np.sum(self.weight_array)
            self.number_of_elements = self.membership.shape[0]
            self.information_content = get_information_content(self.membership)
            if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
    
    cpdef void denoise(self, verbose=False):
        """Examines the membership of each frag. If coverage within a frag
        is < minimum_proportion of coverage around the frag, count the
        overlapping reads as false positive and remove them from the locus."""
        cdef np.ndarray contains_frag, subtable, includes_c, excludes_c, keep
        cdef tuple connections
        cdef Py_ssize_t frag, c_frag
        cdef float includes_weight, excludes_weight, total_weight
        cdef set connected_frags, to_delete
        to_delete = set()
        cdef float removed_weight = 0
        for frag in range(self.membership.shape[1]):
            contains_frag = np.where(self.membership[:,frag]==1)[0]
            if contains_frag.shape[0] > 0:
                # Get a sub-table of reads that include frag
                subtable = self.membership[contains_frag,:]
                subtable[:,frag] = 0
                # Get a table of connected frags that are also members
                connections = np.where(subtable==1)
                connected_frags = set(connections[1])
                for c_frag in connected_frags:
                    # Get all rows where the connected frag was found (1)
                    includes_c = contains_frag[connections[0][connections[1]==c_frag]]
                    # Get all rows where the connected frag is excluded (-1)
                    excludes_c = contains_frag[np.where(subtable[:,c_frag] == -1)[0]]
                    # Compare the weights of the two branches
                    includes_weight = np.sum(self.weight_array[includes_c])
                    excludes_weight = np.sum(self.weight_array[excludes_c])
                    total_weight = includes_weight + excludes_weight
                    if includes_weight < self.minimum_proportion * total_weight: # Excludes strongly outweighs
                        to_delete.update(includes_c)
                        removed_weight += includes_weight
                    elif excludes_weight < self.minimum_proportion * total_weight: # Includes strongly outweighs
                        to_delete.update(excludes_c)
                        removed_weight += excludes_weight

        # Remove all frags that make up < minimum_proportion of their branch
        if verbose:
            if len(to_delete) > 0:
                print('\tDenoising reads: ({} / {}, {}%)'.format(removed_weight, self.weight, round(removed_weight/self.weight*100,1)))
        
        keep = np.array(list(set(range(self.membership.shape[0])).difference(to_delete)))
        if len(keep) > 0:
            self.weight_array = self.weight_array[keep]
            self.rep_array = self.rep_array[keep]
            self.strand_array = self.strand_array[keep]
            self.information_content = self.information_content[keep]
            self.member_content = self.member_content[keep]
            self.membership = self.membership[keep,:]   
            if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
            self.weight = np.sum(self.weight_array)
            self.number_of_elements = self.membership.shape[0]

    cpdef void reduce_membership(self):
        """Given a matrix of membership values, 
        returns a [reduced_membership_matrix, weights] array
        such that all rows are unique and in sort order."""
        cdef np.ndarray reduced_membership, reverse_lookup, new_weights, new_strands, members_bool
        cdef list left_member, right_member, index, sort_triples, sorted_indices
        cdef (int, int, int) triple
        cdef Py_ssize_t i,v
        reduced_membership, reverse_lookup = np.unique(self.membership, axis=0, return_inverse=True)
        new_weights = np.zeros(shape=reduced_membership.shape[0], dtype=np.float32)
        new_strands = np.zeros(shape=reduced_membership.shape[0], dtype=np.int8)
        new_reps = np.zeros(shape=reduced_membership.shape[0], dtype=np.int32)
        
        if type(self) is AnnotationLocus:
            new_traceback = []
            for i in range(reduced_membership.shape[0]):
                new_traceback += [set()]
            
            for i,v in enumerate(reverse_lookup):
                new_traceback[v].add(i)

        for i,v in enumerate(reverse_lookup):
            new_weights[v] += self.weight_array[i]
            new_strands[v] = self.strand_array[i]
            new_reps[v] += self.rep_array[i]

        members_bool = reduced_membership[:,[-4,-1]+list(range(0,reduced_membership.shape[1]-4))+[-3,-2]]==1
        number_of_members = np.sum(members_bool[:,2:-2],axis=1)
        left_member = np.argmax(members_bool, axis=1).tolist()
        right_member = (members_bool.shape[1]-1-np.argmax(members_bool[:,::-1], axis=1)).tolist()
        index = list(range(members_bool.shape[0]))
        sort_triples = sorted(list(zip(left_member, right_member, index)))
        sorted_indices = [triple[2] for triple in sort_triples if number_of_members[triple[2]] > 0]
        self.membership = reduced_membership[sorted_indices,:]
        self.weight_array = new_weights[sorted_indices]
        self.strand_array = new_strands[sorted_indices]
        self.rep_array = new_reps[sorted_indices]
        if len(self.traceback) > 0:
            self.traceback = [new_traceback[i] for i in sorted_indices]

    cpdef void filter_by_reps(self, int minreps=1):
        """Enforce that elements in the membership"""
        cdef np.ndarray keep
        if minreps > 1:
            keep = np.where(self.rep_array >= minreps)[0]
            self.rep_array = self.rep_array[keep]
            self.weight_array = self.weight_array[keep]
            self.strand_array = self.strand_array[keep]
            self.membership = self.membership[keep,:]
            self.number_of_elements = len(keep)
            self.information_content = get_information_content(self.membership)
            if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]

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
            maxIC = self.membership.shape[1]
            new_weights = resolve_containment(self.overlap, self.member_content, self.information_content, self.weight_array, maxIC)
            keep = np.where(new_weights > 0.1)[0]
            self.number_of_elements = len(keep)
            self.overlap = self.overlap[keep,:][:,keep]
            self.membership = self.membership[keep,:]
            self.weight_array = new_weights[keep]
            self.rep_array = self.rep_array[keep]
            self.strand_array = self.strand_array[keep]
            self.information_content = self.information_content[keep]
            self.member_content = self.member_content[keep]
            if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
    
    def build_graph(self, reduce=True):
        """Constructs one or more graphs from 
        connection values (ones) in the overlap matrix.
        Additionally, stores the set of excluded edges for each node as an 'antigraph'
        """
        cdef:
            int i
            char old_strand, new_strand
            list reassigned_coverage
            np.ndarray updated
            float total_coverage
        
        if reduce: # Collapse linear chains prior to graph construction
            updated = np.array([-1])
            while len(updated) > 0:
                updated = self.prune_unreachable_edges()
                self.collapse_linear_chains()
        
        self.graph = ElementGraph(self.overlap, self.membership, self.weight_array, self.strand_array, self.frag_len, self.novelty_ratio)

    cpdef void assemble_transcripts(self, bint complete=False, bint collapse=True):
        self.graph.assemble(self.minimum_proportion, self.weight, self.mean_read_length, self.naive)
        if complete: # Discard all paths that aren't complete
            paths_to_remove = [i for i in range(len(self.graph.paths)) if not self.graph.paths[i].complete]
        else: # Still remove paths with internal gaps
            paths_to_remove = [i for i in range(len(self.graph.paths)) if self.graph.paths[i].has_gaps]
        
        self.graph.remove_paths(sorted(list(set(paths_to_remove))))
        reassigned_coverage = self.reassign_reads_to_paths()
        total_coverage = sum(reassigned_coverage)
        paths_to_remove = self.filter_low_coverage(reassigned_coverage, total_coverage)
        paths_to_remove += self.filter_truncations()
        if self.intron_filter > 0:
            paths_to_remove += self.filter_retained_introns()
        
        self.graph.remove_paths(sorted(list(set(paths_to_remove))))
        reassigned_coverage = self.reassign_reads_to_paths()
        total_coverage = sum(reassigned_coverage)
        paths_to_remove = self.filter_low_coverage(reassigned_coverage, total_coverage)
        while len(paths_to_remove) > 0:
            self.graph.remove_paths(sorted(list(set(paths_to_remove))))
            reassigned_coverage = self.reassign_reads_to_paths()
            total_coverage = sum(reassigned_coverage)
            paths_to_remove = self.filter_low_coverage(reassigned_coverage, total_coverage)
        
        if collapse: # Combine paths whose ends are in range of each other
            pass
            #TODO: Implement finding and combining ends
            # end1.span.stop > end2.span.start

        for i in range(len(self.graph.paths)):
            path = self.graph.paths[i]
            self.transcripts.append(self.convert_path(path, i+1))
        
        self.add_transcript_attributes()
    
    cdef list filter_low_coverage(self, list reassigned_coverage, float total_coverage):
        cdef list low_coverage_paths = [i for i in range(len(reassigned_coverage)) if reassigned_coverage[i] < total_coverage*self.minimum_proportion]
        return low_coverage_paths

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
            junctions = self.BP.J_plus if strand == 1 else self.BP.J_minus
            for junction_hash, junction_count in junctions.items():
                # Convert all intron pair positions to frag indices
                block = self.string_to_span(junction_hash)
                l = block[0] - self.leftmost
                r = block[1] - self.leftmost
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
                                container_reads += other_path.reads
            
            # After all containers are found, compare path.reads to sum of all containers
            if path.reads < container_reads:
                paths_to_remove.append(index)

        return paths_to_remove
    
    cpdef void add_transcript_attributes(self):
        """Populate the new read objects with diagnostic information
        to store in the GTF attributes column."""
        cdef:
            int first, last, s_pos, e_pos
            dict S_info, E_info, SBP, EBP

        for T in self.transcripts:
            T.attributes['length'] = T.get_length()
            T.attributes['reads'] = round(T.weight, 2)
            S_info = {'S.reads':0, 'S.capped':0, 'S.left':0, 'S.right':0}
            E_info = {'E.reads':0, 'E.left':0, 'E.right':0}
            
            if T.strand != 0:
                if T.strand == 1:
                    s_pos = first = T.span[0]
                    e_pos = last = T.span[1]
                    SBP = self.BP.S_plus
                    EBP = self.BP.E_plus
                else:
                    s_pos = first = T.span[1]
                    e_pos = last = T.span[0]
                    SBP = self.BP.S_minus
                    EBP = self.BP.E_minus
                
                if T.s_tag:
                    if first in SBP.keys():
                        bp = self.BP.branchpoints[SBP[first]]
                        S_info['S.reads'] = round(bp.weight,2)
                        S_info['S.capped'] = round(bp.capped,2)
                        S_info['S.left'] = bp.left
                        S_info['S.right'] = bp.right
                    else: # The slow way: Find which branchpoint's span the end is under
                        for bp in self.BP.branchpoints:
                            if bp.branchtype == 'S' and bp.strand == T.strand:
                                if bp.weight > S_info['S.reads']:
                                    if first in bp.span:
                                        S_info['S.reads'] = round(bp.weight,2)
                                        S_info['S.capped'] = round(bp.capped,2)
                                        S_info['S.left'] = bp.left
                                        S_info['S.right'] = bp.right
                                        s_pos = bp.pos

                        if s_pos != first: # S pos was replaced
                            if T.strand == 1:
                                T.ranges[0] = (s_pos, T.ranges[0][1])
                            else:
                                T.ranges[-1] = (T.ranges[-1][0], s_pos)
                
                if T.e_tag:
                    if last in EBP.keys():
                        bp = self.BP.branchpoints[EBP[last]]
                        E_info['E.reads'] = round(bp.weight,2)
                        E_info['E.left'] = bp.left
                        E_info['E.right'] = bp.right
                    else: # The slow way: Find which branchpoint's span the end is under
                        for bp in self.BP.branchpoints:
                            if bp.branchtype == 'E' and bp.strand == T.strand:
                                if bp.weight > E_info['E.reads']:
                                    if last in bp.span:
                                        E_info['E.reads'] = round(bp.weight,2)
                                        E_info['E.left'] = bp.left
                                        E_info['E.right'] = bp.right
                                        e_pos = bp.pos
                    
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
        cdef np.ndarray linear_chains = find_linear_chains(self.overlap)
        cdef dict chain_parent = {}
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
            keep.sort()
            self.number_of_elements = len(keep)
            self.overlap = self.overlap[keep, :][:, keep]
            self.membership = self.membership[keep, :]
            self.strand_array = self.strand_array[keep]
            self.weight_array = self.weight_array[keep]
            self.information_content = self.information_content[keep]
            self.member_content = self.member_content[keep]
            self.reduce_membership()
            self.weight = np.sum(self.weight_array)
            self.number_of_elements = self.membership.shape[0]
            self.information_content = get_information_content(self.membership)
            self.member_content = get_member_content(self.membership)
            self.build_overlap_matrix(True)
            if len(self.traceback) > 0: self.traceback = [self.traceback[k] for k in keep]
    
    cpdef void merge_reads(self, int child_index, int parent_index):
        """Combines the information of two read elements in the locus."""
        cdef char p, c, s
        cdef Py_ssize_t i
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
        
        # Add all of child's weight to parent
        self.weight_array[parent_index] += self.weight_array[child_index]
        self.weight_array[child_index] = 0

    cpdef np.ndarray prune_unreachable_edges(self):
        """Given a overlap matrix, determine which directed edges are unreachable
        from an s-t path on that strand. Removes edges in-place in the matrix.
        """
        cdef int sink_plus, source_minus, source_plus, sink_minus, SOURCE, SINK, a, b, v, w
        cdef char p, m, strand, C
        cdef dict BFS_from_ends, PF, PR, MF, MR
        cdef set from_sink_minus, from_sink_plus, from_source_minus, from_source_plus, in_plus_path, in_minus_path
        cdef np.ndarray new_strands, updated
        cdef char [:, :] COMPATIBILITY = self.overlap
        vertices = COMPATIBILITY.shape[0]
        self.adj = build_adjacencies(self.overlap)
        PF = self.adj['PF']
        PR = self.adj['PR']
        MF = self.adj['MF']
        MR = self.adj['MR']
        new_strands = infer_strands_by_reachability(self.strand_array, self.adj)
        updated = np.where(new_strands != self.strand_array)[0]
        for v in updated: #Check if any strands were updated
            if new_strands[v] == 1: # Strand inferred as plus; remove all minus edges
                self.membership[v,-2:] = -1 
                for w in MF[v]:
                    if new_strands[w] != 0 and new_strands[w] != new_strands[v]:
                        # If the reads are inferred to be on opposite strands, all overlaps are invalidated
                        COMPATIBILITY[v,w] = -1
                        COMPATIBILITY[w,v] = -1
                    else:
                        C = COMPATIBILITY[v, w]
                        if C == 1:
                            COMPATIBILITY[v, w] = 0
                
                for w in MR[v]:
                    if new_strands[w] != 0 and new_strands[w] != new_strands[v]:
                        COMPATIBILITY[v,w] = -1
                        COMPATIBILITY[w,v] = -1
                    else:
                        C = COMPATIBILITY[w, v]
                        if C == 1:
                            COMPATIBILITY[w, v] = 0
            elif new_strands[v] == -1: # Strand inferred as minus; remove all plus edges
                self.membership[v,-4:-2] = -1
                for w in PF[v]:
                    if new_strands[w] != 0 and new_strands[w] != new_strands[v]:
                        # If the reads are inferred to be on opposite strands, all overlaps are invalidated
                        COMPATIBILITY[v,w] = -1
                        COMPATIBILITY[w,v] = -1
                    else:
                        C = COMPATIBILITY[v, w]
                        if C == 1:
                            COMPATIBILITY[v, w] = 0
                
                for w in PR[v]:
                    if new_strands[w] != 0 and new_strands[w] != new_strands[v]:
                        # If the reads are inferred to be on opposite strands, all overlaps are invalidated
                        COMPATIBILITY[v,w] = -1
                        COMPATIBILITY[w,v] = -1
                    else:
                        C = COMPATIBILITY[w, v]
                        if C == 1:
                            COMPATIBILITY[w, v] = 0
        
        self.strand_array = new_strands
        return updated

    cpdef list reassign_reads_to_paths(self):
        """Reassigns reads by setting priors based on unique
        assignments, then distributes shared reads via EM."""
        # Initialize transcripts with uniquely-assignable Elements
        cdef:
            np.ndarray frag_assignment, coverage, priors, priors_count, bin_priors, bin_count
            dict junctions, intron_assignment, intron_assign_count
            int number_of_paths, number_assigned, path_index, i, f, count, uninitialized_count, l, r
            float existing_coverage, f_depth, total_assigned_cov
            (int, int) span, junction
            list frag_assign_count, value, assigned, assigned_cov
        
        coverage = self.BP.depth
        junctions = copy.copy(self.BP.J_plus)
        junctions.update(self.BP.J_minus)
        number_of_paths = len(self.graph.paths)
        number_of_frags = len(self.frags)
        priors = np.zeros(number_of_paths, dtype=np.float32)
        priors_count = np.zeros(number_of_paths, dtype=np.float32)
        
        # Create lookup structures for (partial) exons and introns
        frag_assignment = np.zeros((number_of_paths,number_of_frags), dtype=np.int8)
        intron_assignment = {}
        for i in range(number_of_paths):
            path = self.graph.paths[i]
            for span in path.get_exons():
                frag_assignment[i, span[0]:span[1]+1] = 1
            
            for span in path.get_introns():
                junction = (self.frags[span[0]][0], self.frags[span[1]][1])
                junction_hash = self.span_to_string(junction)
                if junction_hash not in junctions:
                    path.complete = False
                    path.has_gaps = True
                    continue
                
                intron_assignment[junction_hash] = intron_assignment.get(junction_hash, []) + [i]
        
        # Update priors with the coverage depth of all unique features
        frag_assign_count = np.sum(frag_assignment, axis=0).tolist()
        for f in range(number_of_frags):
            count = frag_assign_count[f]
            if count == 1: # Unique assignment
                span = self.frags[f]
                l = span[0] - self.leftmost
                r = span[1] - self.leftmost
                f_depth = np.mean(coverage[l:r])
                for i in np.where(frag_assignment[:,f]==1)[0]:
                    priors[i] += f_depth
                    priors_count[i] += 1
        
        intron_assign_count = {junction_hash:len(value) for junction_hash,value in intron_assignment.items()}
        for junction_hash,count in intron_assign_count.items():
            if count == 1: # Uniquely assigned junction
                f_depth = junctions.get(junction_hash, 0.0)
                for i in intron_assignment[junction_hash]:
                    priors[i] += f_depth
                    priors_count[i] += 1
        
        # Average out the observations of all unique frags/introns
        for i in range(priors.shape[0]):
            if priors_count[i] > 0:
                priors[i] = priors[i]/priors_count[i]
        
        # Update the priors with all frags/introns assigned to >1 path,
        # applying them in bins of increasing count
        bins = sorted(list(set(frag_assign_count)|set(intron_assign_count.values())))
        for b in bins:
            bin_priors = np.zeros(number_of_paths, dtype=np.float32)
            bin_count = np.zeros(number_of_paths, dtype=np.float32)
            for f in range(number_of_frags):
                count = frag_assign_count[f]
                if count == b: # Count matches the current bin
                    span = self.frags[f]
                    l = span[0] - self.leftmost
                    r = span[1] - self.leftmost
                    f_depth = np.mean(coverage[l:r])
                    # Apportion the depth according to the ratio of the priors
                    assigned = np.where(frag_assignment[:,f]==1)[0].tolist()
                    assigned_cov = [priors[i] for i in assigned]
                    uninitialized_count = assigned_cov.count(0)
                    total_assigned_cov = sum(assigned_cov)
                    if total_assigned_cov < f_depth and uninitialized_count > 0:
                        # Uninitialized prior still exists; give remainder
                        remainder = f_depth - total_assigned_cov
                        for i in range(len(assigned_cov)):
                            bin_count[assigned[i]] += 1
                            if assigned_cov[i] == 0:
                                bin_priors[assigned[i]] += remainder/uninitialized_count
                            else:
                                bin_priors[assigned[i]] += assigned_cov[i]
                    elif total_assigned_cov > 0: # Assign all weight according to the ratio of priors
                        for i in range(len(assigned_cov)):
                            bin_priors[assigned[i]] += f_depth*assigned_cov[i]/total_assigned_cov
                            bin_count[assigned[i]] += 1
                    else:
                        for i in range(len(assigned_cov)): # Assign 0 to all
                            bin_count[assigned[i]] += 1

            for junction_hash,count in intron_assign_count.items():
                if count == b: # Count matches the current bin
                    f_depth = junctions.get(junction_hash, 0.0)
                    # Apportion the depth according to the ratio of the priors
                    assigned = intron_assignment[junction_hash]
                    assigned_cov = [priors[i] for i in assigned] 
                    uninitialized_count = assigned_cov.count(0)
                    total_assigned_cov = sum(assigned_cov)
                    if total_assigned_cov < f_depth and uninitialized_count > 0:
                        # Uninitialized prior still exists; give remainder
                        remainder = f_depth - total_assigned_cov
                        for i in range(len(assigned_cov)):
                            bin_count[assigned[i]] += 1
                            if assigned_cov[i] == 0:
                                bin_priors[assigned[i]] += remainder/uninitialized_count
                            else:
                                bin_priors[assigned[i]] += assigned_cov[i]
                    elif total_assigned_cov > 0: # Assign all weight according to the ratio of priors
                        for i in range(len(assigned_cov)):
                            bin_priors[assigned[i]] += f_depth*assigned_cov[i]/total_assigned_cov
                            bin_count[assigned[i]] += 1
                    else:
                        for i in range(len(assigned_cov)): # Assign 0 to all
                            bin_count[assigned[i]] += 1
            
            # Average out the current bin's values with existing priors
            for i in range(priors.shape[0]):
                if bin_count[i] > 0:
                    bin_priors[i] = bin_priors[i]/bin_count[i]
                    priors[i] = (priors[i] + bin_priors[i]) * .5
        
        # After all frags/introns were considered, priors are fully updated.
        # Reassign all elements of the Graph proportionally to the prior ratios.
        for i in range(priors.shape[0]):
            path = self.graph.paths[i]
            path.reads = 0
            path.coverage = priors[i]

        for e in self.graph.elements:
            assigned = e.assigned_to
            assigned_cov = [priors[i] for i in assigned] 
            total_assigned_cov = sum(assigned_cov)
            if total_assigned_cov > 0:
                for i in range(len(assigned_cov)):
                    self.graph.paths[assigned[i]].reads += e.reads*assigned_cov[i]/total_assigned_cov
        
        # return list(priors)
        return list([path.reads for path in self.graph.paths])

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
        gene_id = 'bookend.{}'.format(self.chunk_number)
        transcript_id = 'bookend.{}.{}'.format(self.chunk_number, transcript_number)
        members = sorted(list(element.members))
        nonmembers = element.nonmembers
        N = element.maxIC
        chrom = self.chrom
        source = 0
        strand = element.strand
        ranges = []
        splice = []
        l = r = last_member = -1
        gap_is_splice = True
        junctions = self.BP.J_plus if strand == 1 else self.BP.J_minus
        for m in members: # Convert membership into ranges
            # 0-indexed open doubles
            if m >= len(self.frags):
                break
            
            frag = self.frags[m]
            if last_member == -1: # Uninitialized
                l, r = frag
                # Update leftmost position if it matches an S/E branchpoint
                if element.s_tag and element.strand == 1:
                    if l in self.BP.S_plus:
                        l = self.BP.branchpoints[self.BP.S_plus[l]].pos
                elif element.e_tag and element.strand == -1:
                    if l in self.BP.E_minus:
                        l = self.BP.branchpoints[self.BP.E_minus[l]].pos
            elif last_member == m-1:
                r = frag[1]
            else: # A gap was jumped
                junction_hash = self.span_to_string((r, frag[0]))
                if junction_hash in junctions:
                    gap_is_splice = True
                else:
                    gap_is_splice = False
                
                splice.append(gap_is_splice)
                exon = (l, r)
                ranges.append(exon)
                l, r = frag
            
            last_member = m
        
        # Update rightmost position
        if element.s_tag and element.strand == -1:
            if r in self.BP.S_minus:
                r = self.BP.branchpoints[self.BP.S_minus[r]].pos
        elif element.e_tag and element.strand == 1:
            if r in self.BP.E_plus:
                r = self.BP.branchpoints[self.BP.E_plus[r]].pos

        exon = (l, r)
        ranges.append(exon)
        s_tag = element.s_tag
        e_tag = element.e_tag
        capped = False
        weight = element.reads
        elementAttributes = {}
        elementData = ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)
        readObject = RNAseqMapping(elementData, elementAttributes)
        readObject.attributes['gene_id'] = gene_id
        readObject.attributes['transcript_id'] = transcript_id
        if element.coverage == -1:
            readObject.coverage = element.mean_coverage(self.mean_read_length)
        else:
            readObject.coverage = element.coverage
        
        readObject.attributes['cov'] = round(readObject.coverage, 2)
        return readObject
    
    def dump_json(self):
        """Returns a string in JSON format that fully describes the locus,
        its input data and its solution."""
        locusname = 'bookend.{}'.format(self.chunk_number)

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
        self.minreps = minreps
        self.confidence = confidence
        self.ignore_reference_ends = ignore_reference_ends
        self.cap_percent = cap_percent
        self.ref_reads, nonref_reads = self.get_annotation_info(list_of_reads)
        if len(nonref_reads) > 0:
            Locus.__init__(self, chrom, chunk_number, nonref_reads, 0, end_extend, min_overhang, True, minimum_proportion, 0, 1, False, False, True, intron_filter, False, False, True)
            self.cap_percent = cap_percent
            self.total_s = sum([bp.weight for bp in self.BP.branchpoints if bp.branchtype == 'S'])
            self.total_e = sum([bp.weight for bp in self.BP.branchpoints if bp.branchtype == 'E'])
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
        self.weight_array = self.weight_array[keep]
        self.rep_array = self.rep_array[keep]
        self.strand_array = self.strand_array[keep]
        self.information_content = self.information_content[keep]
        self.member_content = self.member_content[keep]
        self.transcripts = [self.transcripts[k] for k in keep]
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
        elementData = ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)
        readObject = RNAseqMapping(elementData, elementAttributes)
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

cpdef np.ndarray resolve_containment(np.ndarray overlap_matrix, np.ndarray member_content, np.ndarray information_content, np.ndarray read_weights, int maxIC):
    """Given a overlap matrix, 'bubble up' the weight of
    all reads that have one or more 'contained by' relationships
    to other reads. Pass from highest complexity reads down, assigning
    weight proportional to the existing weight.
    The resulting matrix should contain only overlaps, exclusions, and unknowns."""
    cdef:
        np.ndarray containment, contained, IC_order, new_weights, containers, incompatible_with_containers, incompatible_with_i, incompatible
        Py_ssize_t full_path, i
        float incompatible_weight, total
    
    containment = overlap_matrix==2 # Make a boolean matrix of which reads are contained in other reads
    np.put(containment, range(0,containment.shape[0]**2,containment.shape[0]+1), False, mode='wrap') # Blank out the diagonal (self-containments)
    for full_path in np.where(information_content == maxIC)[0]: # Ignore full paths here
        containment[:,full_path] = False
    
    contained = np.where(np.sum(containment, axis=1) > 0)[0] # Identify reads that are contained
    IC_order = contained[np.lexsort((-information_content[contained], -member_content[contained]))] # Rank them by decreasing number of members
    new_weights = np.copy(read_weights)
    for i in IC_order:
        containers = np.where(containment[i,:])[0]
        
        # Get the set of reads incompatible with all containers but that do not exclude i
        incompatible_with_containers = overlap_matrix[:,containers]==-1
        incompatible_with_i = np.where(overlap_matrix[i,:]==-1)[0]
        incompatible_with_containers[incompatible_with_i,:] = False
        incompatible = np.where(np.all(incompatible_with_containers,axis=1))[0]

        # Calculate the proportion of i to merge into i's containers and the proportion to separate
        if len(incompatible) > 0:
            incompatible_weight = np.sum(new_weights[incompatible])
            # print("[{}] {} is incompatible with containers {}".format(i, incompatible, containers))
        else:
            incompatible_weight = 0
        
        adjusted_weights = new_weights[containers]
        total = np.sum(adjusted_weights) + incompatible_weight
        new_weights[containers] += (new_weights[i] / total) * adjusted_weights
        if incompatible_weight == 0:
            new_weights[i] = 0
            containment[:,i] = False
        else:
            new_weights[i] = incompatible_weight/total * new_weights[i]

    return new_weights


cpdef set dictBFS(dict adjacency, source):
    """Given an adjacency dict, perform Breadth-First Search and return a list of keys that were visited"""
    visited = { v:False for v in adjacency.keys() }
    visited[source] = True
    queue = deque(maxlen=len(visited))
    queue.append(source)
    while queue:
        v = queue.popleft()
        for w in adjacency[v]:
            if not visited[w]:
                visited[w] = True
                queue.append(w)
    
    return set([v for v in visited.keys() if visited[v]])

cpdef dict build_adjacencies(np.ndarray[char, ndim=2] overlap):
    """Returns a dictionary of four adjacency-list graphs,
    a forward and reverse for each strand."""
    cdef Py_ssize_t vertices, a, b
    cdef char p, m
    cdef char [:, :] COMPATIBILITY = overlap
    vertices = COMPATIBILITY.shape[0]
    cdef dict plus_forward = { v:[] for v in range(vertices+4) }
    cdef dict minus_forward = { v:[] for v in range(vertices+4) }
    cdef dict plus_reverse = { v:[] for v in range(vertices+4) }
    cdef dict minus_reverse = { v:[] for v in range(vertices+4) }
    for a in range(vertices-1):
        for b in range(a+1, vertices):
            p = COMPATIBILITY[a, b] # plus edge
            m = COMPATIBILITY[b, a] # minus edge
            if p >= 1:
                plus_forward[a].append(b)
                plus_reverse[b].append(a)
            
            if m >= 1:
                minus_forward[b].append(a)
                minus_reverse[a].append(b)

    return {'PF':plus_forward, 'PR':plus_reverse, 'MF':minus_forward, 'MR':minus_reverse}

cpdef np.ndarray infer_strands_by_reachability(np.ndarray[char, ndim=1] strand_array, dict adj):
    """Given a set of adjacency matrices from prune_unreachable_edges(),
    assigns strands to all nonstranded elements that are only reachable
    from unambiguously stranded reads of a single direction."""
    cdef set reachable_from_plus, reachable_from_minus
    cdef int PSN, MSN, node, min_plus, min_minus, max_plus, max_minus, r
    cdef np.ndarray plus_stranded_nodes = np.where(strand_array == 1)[0]
    cdef np.ndarray minus_stranded_nodes = np.where(strand_array == -1)[0]
    cdef np.ndarray output_strand_array = np.copy(strand_array)
    # Special case: Strands beyond the boundaries of the stranded elements
    # must be reachable from the last stranded element
    cdef dict PF, PR, MF, MR
    PF = adj['PF']
    PR = adj['PR']
    MF = adj['MF']
    MR = adj['MR']
    reachable_from_plus = set()
    reachable_from_minus = set()
    if len(plus_stranded_nodes) > 0:
        min_plus = min(plus_stranded_nodes)
        max_plus = max(plus_stranded_nodes)
        left_of_plus = dictBFS(PR, min_plus)
        right_of_plus = dictBFS(PF, max_plus)
        for PSN in plus_stranded_nodes:
            if PSN not in reachable_from_plus:
                if PSN > min_plus and PSN < max_plus:
                    reachable_from_plus.update(dictBFS(PF, PSN) | dictBFS(PR, PSN))
        
        reachable_from_plus = set([r for r in reachable_from_plus if (r > min_plus or r in left_of_plus) and (r < max_plus or r in right_of_plus)])
    
    if len(minus_stranded_nodes) > 0:
        min_minus = min(minus_stranded_nodes)
        max_minus = max(minus_stranded_nodes)
        left_of_minus = dictBFS(MF, min_minus)
        right_of_minus = dictBFS(MR, max_minus)
        for MSN in minus_stranded_nodes:
            if MSN not in reachable_from_minus:
                if MSN > min_minus and MSN < max_minus:
                    reachable_from_minus.update(dictBFS(MF, MSN) | dictBFS(MR, MSN))
        
        reachable_from_minus = set([r for r in reachable_from_minus if (r > min_minus or r in left_of_minus) and (r < max_minus or r in right_of_minus)])
    
    for node in reachable_from_plus.difference(reachable_from_minus):
        output_strand_array[node] = 1
    
    for node in reachable_from_minus.difference(reachable_from_plus):
        output_strand_array[node] = -1
    
    return output_strand_array

cpdef np.ndarray find_linear_chains(np.ndarray[char, ndim=2] overlap_matrix):
    cdef:
        np.ndarray edges, ingroup, outgroup, in_chain, putative_chain_starts
        int chain, v, next_v
    
    edges = overlap_matrix >= 1
    np.put(edges, range(0,edges.shape[0]**2,edges.shape[0]+1), False, mode='wrap') # Blank out the diagonal (self-containments)
    ingroup = np.sum(edges, axis=0, dtype=np.int32)
    outgroup = np.sum(edges, axis=1, dtype=np.int32)
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
                if outgroup[next_v] != 1:
                    break
                
                next_v = np.where(edges[next_v,:])[0][0]
                if in_chain[next_v] == chain: # Already traversed
                    break
    
    return in_chain


def sum_subset(mask, array_to_mask):
    """Given an array and a boolean mask, return
    the sum of array values at which mask was True"""
    return array_to_mask[mask].sum()

