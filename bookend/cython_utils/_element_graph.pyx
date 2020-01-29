#cython: language_level=3
import json
import cython
import numpy as np
cimport numpy as np
import copy
from networkx import DiGraph
from collections import deque

inf = float('Inf')

cdef class ElementGraph:
    cdef public list elements, paths
    cdef public np.ndarray assignments
    cdef public float novelty_ratio
    cdef readonly int number_of_elements
    cdef Element emptyPath
    def __init__(self, np.ndarray overlap_matrix, np.ndarray membership_matrix, weights, strands, lengths, float novelty_ratio=1):
        """Constructs a forward and reverse directed graph from the
        connection values (ones) in the overlap matrix.
        Additionally, stores the set of excluded edges for each node as an 'antigraph'
        """
        self.emptyPath = Element(-1, 0, 0, np.array([0]), np.array([0]), np.array([0]))
        cdef Element e, path, part
        cdef int e_index, path_index, i
        self.novelty_ratio = novelty_ratio
        self.number_of_elements = overlap_matrix.shape[0]
        self.elements = [Element(
            i, weights[i], strands[i],
            membership_matrix, overlap_matrix, lengths
        ) for i in range(self.number_of_elements)] # Generate an array of Element objects
        self.assignments = np.zeros(shape=self.number_of_elements, dtype=np.int32)
        self.paths = []
        # Assign all reads to any existing complete paths
        for e in self.elements:
            if e.complete: # A full-length path exists in the input elements
                path = copy.deepcopy(e)
                path_index = len(self.paths)
                self.paths.append(path)
                contained = np.where(overlap_matrix[:,e.index]==2)[0]
                for i in contained:
                    if i != path.index:
                        path.includes.add(i)
                        part = self.elements[i]
                        path.reads += self.available_reads(path, part)
                        part.assigned_to.append(path_index)
                        self.assignments[i] += 1
                        part.update()
                    else: # Only update assignments
                        self.assignments[i] += 1
                        part = self.elements[i]
                        part.assigned_to.append(path_index)
                        part.update()
                
                path.update()

    cpdef void assemble(self, float minimum_proportion, float total_reads, float mean_read_length, bint naive):
        """Iteratively perform find_optimal_path() on the graph
        until the number of novel reads fails to exceed minimum_proportion
        of the reads at the locus. If minimum_proportion == 0, assemble()
        only terminates when every read is in a path."""
        cdef float threshold, total_reads_assigned
        cdef Element path
        
        total_reads_assigned = sum([path.reads for path in self.paths])
        path = self.find_optimal_path(mean_read_length, naive)
        if path is self.emptyPath:
            return
        
        threshold = total_reads*(1-minimum_proportion)
        total_reads_assigned += self.add_path(path)
        while total_reads_assigned < threshold:
            path = self.find_optimal_path(mean_read_length, naive)
            if path is self.emptyPath:
                total_reads_assigned = threshold
            else:
                total_reads_assigned += self.add_path(path)

    cpdef float available_reads(self, Element path, Element part):
        """Given a path that wants to merge with the indexed element,
        calculate how much coverage is actually available to the path."""
        # Get the total cov of all already assigned paths
        if len(part.assigned_to) == 0: # No competition, all reads are available
            return part.reads
        
        cdef float assigned_cov = path.cov
        for i in part.assigned_to:
            assigned_cov += self.paths[i].cov
        
        proportion = path.cov / assigned_cov
        return proportion * part.reads
    
    cpdef int max_info_gained(self, Element fromElement, Element throughElement):
        """Returns the amount of unique information in the throughElement."""
        new_members = len(throughElement.members.difference(fromElement.members))
        new_nonmembers = len(throughElement.nonmembers.difference(fromElement.nonmembers))
        
        return new_members+new_nonmembers

    cpdef void take_from_competitors(self, Element e, float proportion, bint naive):
        """Remove the proportion of weight from competing paths."""
        cdef float total_cov, c, partial_proportion, reads_to_remove
        cdef int i, num_competitors
        cdef list competitor_cov, competitor_proportion
        num_competitors = len(e.assigned_to)
        if num_competitors > 0:
            if naive:
                # NAIVE: Remove an equal amount or reads from each competitor
                partial_proportion = proportion*e.reads/num_competitors
                for i in range(num_competitors):
                    self.paths[e.assigned_to[i]].reads -= partial_proportion
            else:
                # PROPORTIONAL: Each competitor is "selfish" and holds onto a proportional number of reads
                competitor_cov = [self.elements[i].cov for i in e.assigned_to]
                total_cov = sum(competitor_cov)
                competitor_proportion = [c/total_cov for c in competitor_cov]
                reads_to_remove = e.reads*proportion
                for i in range(num_competitors):
                    self.paths[e.assigned_to[i]].reads -= reads_to_remove*competitor_proportion[i]
    
    cpdef void add_edge_to_path(self, Element path, dict edge, bint naive):
        """Merges the proper 
        """
        cdef float proportion = edge['available'] / edge['total']
        path.merge(edge['element'], proportion)
        self.take_from_competitors(edge['element'], proportion, naive)  

    cpdef list pick_best_edge(self, Element path, dict edge_dict, float mean_read_length):
        """When one or more mutually incompatible edges are available to traverse,
        the path branches and a decision must be made about which branch to follow.
        For each edge, a score is calculated:
            
        """
        cdef Element element, other_element, e
        cdef set compatible
        cdef int c, i, j, length_gained, info_gained, min_gained
        cdef list best_indices, extend_list, included_elements
        cdef dict max_extend
        cdef float reads_gained, score, best_score, dead_end_penalty
        cdef bint is_internal
        dead_end_penalty = 0.1
        best_indices = []
        best_score = 0
        min_gained = max(int(mean_read_length*.5),1)
        # Choose edges by decreasing absolute distance from path.index
        elements = [edge_dict[i]['element'] for i in edge_dict.keys()]
        #TODO: Fix decision tree for when elements are INSIDE path (filling a gap)
        extend_left = {}
        extend_right = {}
        max_extend = {}
        for e in elements:
            extend_left[e.index] = path.LM - e.LM 
            extend_right[e.index] = e.RM - path.RM
            max_extend[e.index] = max([extend_left[e.index], extend_right[e.index]])
        
        extend_list = sorted([(v,k) for k,v in max_extend.items()], reverse=True)
        for m,i in extend_list:
            is_internal = m < 0
            element = edge_dict[i]['element']
            # Build a list of edge indices reciprocally compatible with edge
            compatible = set([i])
            compatible.update([ # Populate the set of compatible elements with those that...
                e.index for e in elements # Are in the list of elements,
                if e.strand == element.strand # Are on the same strand,
                and max_extend[e.index] > 0 == is_internal # Are both inside or both outside of path,
                and extend_left[e.index] < extend_left[element.index]  # Doesn't extend past element on the left
                and extend_right[e.index] < extend_right[element.index] # Doesn't extend past element on the right
            ])
            compatible.difference_update(element.excludes)
            included_elements = sorted([self.elements[c] for c in compatible], reverse=True)
            for e in included_elements: # From heaviest to lightest, exclude all incompatibilities of included elements
                if e.index in compatible:
                    compatible.difference_update(e.excludes)
            
            reads_gained = sum([edge_dict[c]['available'] for c in compatible])
            length_gained = max(edge_dict[i]['length'], min_gained)
            info_gained = 1
            # info_gained = edge_dict[i]['newinfo']
            # if element.is_spliced and mean_read_length < length_gained:
                # length_gained = round(mean_read_length)
            
            # score = info_gained * reads_gained / ((length_gained+1) * ((len(element.assigned_to)*self.novelty_ratio)+1))
            score = reads_gained / length_gained
            if self.is_dead_end(path, element):
                score = score * dead_end_penalty
            
            if score > best_score:
                best_score = score
                best_indices = list(compatible)
        
        return best_indices

    cpdef Element get_heaviest_element(self):
        cdef Element best_element, new_element
        cdef np.ndarray available_elements = np.where(self.assignments==0)[0]
        if len(available_elements) == 0:
            return self.emptyPath
        
        best_element = self.elements[available_elements[0]]
        for i in available_elements:
            new_element = self.elements[i]
            if new_element > best_element:
                best_element = new_element
            elif new_element == best_element: # Break ties by complexity
                if new_element.IC > best_element.IC:
                    best_element = new_element
        
        return copy.deepcopy(best_element)

    cpdef dict make_edge(self, Element path, Element e):
        return {
            'element':e,
            'total':e.reads, 
            'available':self.available_reads(path, e), 
            'length':e.uniqueLength(path),
            'newinfo':e.uniqueInformation(path)
        }
    
    cpdef bint is_dead_end(self, Element path, Element next_path):
        """Returns True iff path merging with next_path would
        result in no available edges in the direction(s) being extended
        and next_path is not a valid terminator in that direction."""
        cdef Element downstream_path
        cdef list available_from_next_path
        cdef bint extending_left, extending_right
        if next_path.s_tag or next_path.e_tag: # True ends cannot be dead ends
            return False
        
        extending_left = extending_right = False
        if next_path.LM < path.LM:
            extending_left = True
        
        if next_path.RM > path.RM:
            extending_right = True
        
        #TODO: Redefine dead ends so terminating left OR right counts as a dead end when extending in two directions
        available_from_next_path = [self.elements[i] for i in next_path.ingroup | next_path.outgroup if i != path.index]
        for downstream_path in available_from_next_path:
            if downstream_path.compatible(path): # At least one path connected to next_path can continue path
                if downstream_path.s_tag or downstream_path.e_tag: # Path can continue to a terminal element
                    return False
                
                if extending_left and downstream_path.LM < next_path.LM: # Path can continue left
                    return False
                
                if extending_right and downstream_path.RM > next_path.RM: # Path can continue right
                    return False
        
        return True # At least one extension is still possible if path and next_path merge
            

    cpdef Element find_optimal_path(self, float mean_read_length, bint naive, bint backtrack=True, bint step=False):
        """Traverses the path in a greedy fashion from the heaviest element."""
        cdef Element currentPath, e
        cdef int i, counter
        cdef list edges, best_indices
        cdef set available
        cdef dict edge_dict, edge
        cdef np.ndarray available_elements
        # Get the current working path (heaviest unassigned Element)
        currentPath = self.get_heaviest_element()
        if currentPath is self.emptyPath:
            return currentPath
        
        available = currentPath.ingroup | currentPath.outgroup
        # while len(available) > 0: # Extend as long as possible
        for counter in range(1000):
            if len(available) == 0:
                break
            
            best_indices = []
            if len(available) == 1: # Only one option, do not evaluate
                i = list(available)[0]
                e = self.elements[i]
                edge = self.make_edge(currentPath, e)
                self.add_edge_to_path(currentPath, edge, naive)
            else: # Pick the heaviest edge
                edges = sorted(list(available))
                edge_dict = {i:self.make_edge(currentPath, self.elements[i]) for i in edges}
                if all_are_compatible([edge['element'] for edge in edge_dict.values()]):
                    best_indices = edges
                else:
                    free_edges = [i for i in edge_dict.keys() if edge_dict[i]['newinfo'] == 0]
                    for i in free_edges: # Merge all edges that don't change the structure
                        e = edge_dict[i]['element']
                        if e.strand != 0 and currentPath.strand != e.strand: # Not actually free, would add strandedness
                            continue
                        
                        self.add_edge_to_path(currentPath, edge_dict[i], naive)
                        del edge_dict[i]
                    
                    if len(edge_dict) == 1:
                        edge = list(edge_dict.values())[0]
                        self.add_edge_to_path(currentPath, edge, naive)
                    elif len(edge_dict) > 0:
                        best_indices = self.pick_best_edge(currentPath, edge_dict, mean_read_length)
                        if len(best_indices) == 0:
                            print("ERROR, no edge chosen!")
                
                if step:
                    print(currentPath)
                    print("Merging elements {}:".format(best_indices))
                    for i in best_indices:
                        print(self.elements[i])
                    
                    input('...')
                
                for i in best_indices:
                    # Add the element to the path with a proportional amount of reads
                    if i in edge_dict:
                        edge = edge_dict[i]
                    else: # A compatibility was found that wasn't in the original set
                        e = self.elements[i]
                        edge = self.make_edge(currentPath, e)
                    
                    self.add_edge_to_path(currentPath, edge, naive)
            
            available = currentPath.ingroup | currentPath.outgroup
        
        if counter == 999:
            print("ERROR: stuck in loop. Breaking.")
        
        # A terminating path was found. Add all compatible reads that were missed
        for i in range(len(self.elements)):
            e = self.elements[i]
            if currentPath.LM <= e.LM and currentPath.RM >= e.RM: # E is contained in the bounds of currentPath
                if e.index not in currentPath.excludes: # E not excluded by currentPath
                    if e.index not in currentPath.includes: # E not included in currentPath
                        if e.compatible(currentPath):
                            edge = {
                                'element':e,
                                'total':e.reads, 
                                'available':self.available_reads(currentPath, e), 
                                'length':e.uniqueLength(currentPath)
                            }
                            self.add_edge_to_path(currentPath, edge, naive)
        
        return currentPath

    cpdef float add_path(self, Element path):
        """Evaluate what proportion of the compatible reads should be """
        cdef Element e, p
        cdef int i
        cdef float competing_reads, percent_in_path, extraction_depth, extracted_reads, portion_to_extract, amount_to_extract
        cdef set other_edges = set()
        cdef float novel_reads = 0
        # Assign each included element to the path
        for i in path.includes:
            if self.assignments[i] == 0:
                novel_reads += self.elements[i].reads
            
            self.assignments[i] += 1
            self.elements[i].assigned_to.append(len(self.paths))
            self.elements[i].update()
        
        # Add the new path to the list of paths
        self.paths.append(path)
        return novel_reads
    
    cpdef void remove_paths(self, list indices):
        """Removes all trace of a path from paths."""
        cdef Element path, element
        cdef int i, index
        cdef list keep
        cdef dict old_indices = {i:self.paths[i] for i in range(len(self.paths))}
        if len(indices) == 0:
            return
        
        for index in indices:
            path = self.paths[index]
            for i in path.includes:
                element = self.elements[i]
                element.assigned_to.remove(index)
                self.assignments[i]-=1
        
        keep = [index for index in range(len(self.paths)) if index not in indices]
        self.paths = [self.paths[i] for i in keep]
        for i in range(len(self.paths)): # Update the index attribute of each path
            self.paths[i].index = i
        
        for i in range(len(self.elements)): # Update each assigned_to to keep elements connected to paths
            element = self.elements[i]
            element.assigned_to = [old_indices[a].index for a in element.assigned_to]


###########################################

cdef class Element:
    """Represents a read or collection of reads in a Locus."""
    cdef public int index, length, IC, maxIC, left, right, number_of_elements, LM, RM
    cdef public list assigned_to
    cdef public char strand
    cdef public float cov, reads, coverage
    cdef public set members, nonmembers, ingroup, outgroup, excludes, includes, end_indices
    cdef public np.ndarray frag_len
    cdef public bint complete, s_tag, e_tag, empty, is_spliced, has_gaps
    def __init__(self, int index, float reads, char strand, np.ndarray membership, np.ndarray overlap, np.ndarray frag_len):
        cdef Py_ssize_t i
        cdef char m, overOut, overIn
        self.is_spliced = False
        self.index = self.left = self.right = index
        self.frag_len = frag_len
        self.number_of_elements = overlap.shape[0]
        self.includes = set([self.index])
        self.reads = reads
        self.strand = strand
        self.length = 0
        self.coverage = -1
        self.assigned_to = []
        self.complete = False
        self.has_gaps = False
        self.members = set()
        self.nonmembers = set()
        self.ingroup = set()
        self.outgroup = set()
        self.excludes = set()
        if index == -1: # Special Element emptyPath: placeholder for null values
            self.empty = True
            self.maxIC = 0
            self.end_indices = set()
        else:
            self.empty = False
            self.maxIC = membership.shape[1]
            self.end_indices = set(range(self.maxIC-4, self.maxIC))
            for i in range(membership.shape[1]):
                m = membership[self.index, i]
                if m == 1:
                    self.members.add(i)
                    self.length += self.frag_len[i]
                elif m == -1:
                    self.nonmembers.add(i)
                else:
                    continue
            
            self.update()
            for i in range(self.number_of_elements):
                if i == self.index:
                    continue
                
                overOut = overlap[self.index, i]
                if overOut == -1:
                    self.excludes.add(i)
                    continue
                elif overOut >= 1:
                    self.outgroup.add(i)
                
                overIn = overlap[i, self.index]
                if overIn >= 1:
                    self.ingroup.add(i)
    
    def __repr__(self):
        chars = [' ']*self.maxIC
        strand = {-1:'-', 0:'.', 1:'+'}[self.strand]
        for m in self.members:
            chars[m] = '*'
        
        for n in self.nonmembers:
            chars[n] = '_'
        
        return '|{}| {}-{} ({})'.format(''.join(chars),self.left,self.right,strand)

    cpdef float mean_coverage(self, float mean_read_length):
        """Returns an estimated coverage depth per nucleotide."""
        cdef float estimated_coverage = self.reads * mean_read_length
        return estimated_coverage / self.length

    cpdef str as_string(self):
        cdef str string = ''
        cdef int i
        for i in range(self.number_of_elements):
            if i in self.includes:
                string+='+'
            elif i in self.excludes:
                string+='-'
            else:
                string+=' '
        
        return string
    
    def __eq__(self, other): return self.cov == other.cov
    def __ne__(self, other): return self.cov != other.cov
    def __gt__(self, other): return self.cov >  other.cov
    def __ge__(self, other): return self.cov >= other.cov
    def __lt__(self, other): return self.cov <  other.cov
    def __le__(self, other): return self.cov <= other.cov

    def __add__(self, other):
        if self.empty: # The special cast emptyPath defeats addition
            return self
        elif other.empty:
            return other
        
        summed_element = copy.deepcopy(self)
        if other.index in self.outgroup:
            forward = True
        elif other.index in self.ingroup:
            forward = False
        else:
            raise Exception('Error: Element {} is not connected to Element {}'.format(other, self))
        
        summed_element.merge(other, 1)
        return summed_element
    
    cpdef void update(self):
        cdef int n
        if self.empty:
            return
        
        self.LM = min(self.members.difference(self.end_indices))
        self.RM = max(self.members.difference(self.end_indices))
        if self.strand == 1:
            self.s_tag = self.maxIC - 4 in self.members # has + start
            self.e_tag = self.maxIC - 3 in self.members # has + end
        elif self.strand == -1:
            self.s_tag = self.maxIC - 2 in self.members # has + start
            self.e_tag = self.maxIC - 1 in self.members # has + end
        
        self.cov = self.reads / (self.length * (len(self.assigned_to)+1))
        self.IC = len(self.members) + len(self.nonmembers)
        if self.IC == self.maxIC:
            self.complete = True
        else:
            self.complete = False
        
        for n in self.nonmembers:
            if n > self.LM and n < self.RM:
                self.is_spliced = True
                break
    
    cpdef set uniqueMembers(self, Element other):
        """Given a second Element, return a set of frags that
        are only in self and not in other"""
        return self.members.difference(other.members)
    
    cpdef set uniqueNonembers(self, Element other):
        """Given a second Element, return a set of frags that
        are only in self and not in other"""
        return self.nonmembers.difference(other.nonmembers)

    cpdef int uniqueInformation(self, Element other):
        """Given a second Element, return a set of frags that
        are only in self and not in other"""
        return len(self.uniqueMembers(other)|self.uniqueNonembers(other))

    cpdef int uniqueLength(self, Element other):
        """Given a second Element, return the total length that is unique
        to self (summed length of uniqueMembers)."""
        cdef int length = 0
        for m in self.uniqueMembers(other):
            length += self.frag_len[m]
        
        return length

    cpdef bint compatible(self, Element other):
        """Returns a boolean of whether or not self and other could be
        subpaths in a shared path."""
        if self.empty: # emptyPath is incompatible with everything
            return False
        
        if self.strand != 0 and other.strand != 0 and self.strand != other.strand:
            # self and other are on opposite strands
            return False
        
        if self.members.isdisjoint(other.nonmembers): # Must not contain any excluded frags
            if other.members.isdisjoint(self.nonmembers): # Check reciprocal
                return True
        
        return False
    
    cpdef list get_introns(self):
        """Returns a list of membership index pairs that mark the
        start and end of each internal gap (intron) in the path."""
        cdef:
            int istart, iend, i, n
            list introns, internal_skips
        
        introns = []
        internal_skips = sorted([n for n in self.nonmembers if n > self.LM and n < self.RM])
        istart = iend = -1
        for i in internal_skips:
            if iend == -1: # Uninitialized
                istart = iend = i
            elif i == iend + 1: # Contiguous
                iend = i
            else: # Gapped
                introns.append((istart,iend))
                istart = iend = i
        
        if istart != -1:
            introns.append((istart,iend))
        
        return introns
    
    cpdef list get_exons(self):
        """Returns a list of membership index pairs that mark the
        start and end of each contiguous stretch (exon) in the path."""
        cdef:
            int istart, iend, i, m
            list exons, internal_members
        
        exons = []
        internal_members = sorted([m for m in self.members if m >= self.LM and m <= self.RM])
        istart = iend = -1
        for i in internal_members:
            if iend == -1: # Uninitialized
                istart = iend = i
            elif i == iend + 1: # Contiguous
                iend = i
            else: # Gapped
                exons.append((istart,iend))
                istart = iend = i
        
        if istart != -1:
            exons.append((istart,iend))
        
        return exons

    cpdef void merge(self, Element other, float proportion):
        """Add an Element to this one, combining their membership and reads
        through in-place updates of self."""
        cdef set covered, unique
        cdef bint extendedLeft, extendedRight
        cdef int o, i, leftMember, rightMember
        if self.empty:
            return
        
        if not self.compatible(other):
            print('ERROR: {} incompatible with {}'.format(self, other))
            self.outgroup.discard(other.index)
            self.ingroup.discard(other.index)
            other.outgroup.discard(self.index)
            other.ingroup.discard(self.index)
            return
        
        if self.strand == 0:
            self.strand = other.strand
        
        self.reads += other.reads * proportion
        unique = other.uniqueMembers(self)
        # Update Membership
        self.nonmembers.update(other.nonmembers)
        # Update Overlaps
        self.excludes.update(other.excludes) # Sum of exclusions
        self.includes.update(other.includes) # Sum of inclusions
        self.outgroup.difference_update(self.excludes)
        self.ingroup.difference_update(self.excludes)
        if len(unique) > 0: 
            leftMember = min(self.members.difference(self.end_indices))
            rightMember = max(self.members.difference(self.end_indices))
            extendedLeft = False
            extendedRight = False
            for f in unique: # Append each frag from other to self
                self.members.add(f)
                self.length += self.frag_len[f]
                if f not in self.end_indices:
                    if f > rightMember: 
                        extendedRight = True
                    elif f < leftMember:
                        extendedLeft = True
            
            if extendedRight: # Replace the rightward edges with those of Other
                self.outgroup = set([o for o in self.outgroup if o < self.left])
                self.ingroup = set([i for i in self.ingroup if i < self.left])
            
            if extendedLeft: # Replace the leftward edges with those of Other
                self.outgroup = set([o for o in self.outgroup if o > self.right])
                self.ingroup = set([i for i in self.ingroup if i > self.right])
            
            self.outgroup.update(other.outgroup.difference(self.excludes))
            self.ingroup.update(other.ingroup.difference(self.excludes))
            
            # Update the left and right borders of the Element
            self.right = max(self.right, other.right)
            self.left = min(self.left, other.left)
            covered = set(range(self.left,self.right))
            self.outgroup.difference_update(covered)
            self.ingroup.difference_update(covered)
        else:
            # Update the left and right borders of the Element
            self.right = max(self.right, other.right)
            self.left = min(self.left, other.left)
        
        self.update()
        if self.strand == 1: # Enforce directionality of edges
            self.outgroup = set([o for o in self.outgroup if o > self.right])
            self.ingroup = set([i for i in self.ingroup if i < self.left])
        elif self.strand == -1:
            self.outgroup = set([o for o in self.outgroup if o < self.left])
            self.ingroup = set([i for i in self.ingroup if i > self.right])

#########################################            

cpdef bint all_are_compatible(list element_list):
    """Iterates over a list of elements, returning whether
    or not all Elements are mutually compatible."""
    cdef Element e
    cdef set incompatible
    cdef set excludes = set()
    for e in element_list:
        incompatible = e.includes.intersection(excludes)
        if len(incompatible) > 0:
            return False
        
        excludes.update(e.excludes)
    
    return True
