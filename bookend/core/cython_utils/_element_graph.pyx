#cython: language_level=3
import json
import cython
import numpy as np
cimport numpy as np
import copy
#from networkx import DiGraph
from collections import deque

inf = float('Inf')

cdef class ElementGraph:
    cdef public list elements, paths
    cdef public np.ndarray assignments, overlap
    cdef readonly int number_of_elements, maxIC
    cdef Element emptyPath
    cdef public float bases
    cdef public set SP, SM, EP, EM
    cdef public bint no_ends
    def __init__(self, np.ndarray overlap_matrix, np.ndarray membership_matrix, weights, strands, lengths):
        """Constructs a forward and reverse directed graph from the
        connection values (ones) in the overlap matrix.
        Additionally, stores the set of excluded edges for each node as an 'antigraph'
        """
        self.emptyPath = Element(-1, np.array([0]), 0, np.array([0]), np.array([0]), np.array([0]), 0)
        cdef Element e, path, part
        cdef int e_index, path_index, i
        self.SP, self.SM, self.EP, self.EM = set(), set(), set(), set()
        self.overlap = overlap_matrix
        self.number_of_elements = self.overlap.shape[0]
        self.maxIC = membership_matrix.shape[1]
        self.elements = [Element(
            i, weights[i,:], strands[i],
            membership_matrix[i,:], self.overlap, lengths, self.maxIC
        ) for i in range(self.number_of_elements)] # Generate an array of Element objects
        self.assignments = np.zeros(shape=self.number_of_elements, dtype=np.int32)
        self.paths = []
        self.bases = sum([e.bases for e in self.elements])
        self.check_for_full_paths()
    
    cpdef void check_for_full_paths(self):
        """Assign all reads to any existing complete paths"""
        cdef Element e, path, part
        cdef np.ndarray contained
        for e in self.elements:
            if e.strand == 1:
                if e.s_tag: self.SP.add(e.index)
                if e.e_tag: self.EP.add(e.index)
            elif e.strand == -1:
                if e.s_tag: self.SM.add(e.index)
                if e.e_tag: self.EM.add(e.index)
            if e.complete: # A full-length path exists in the input elements
                path = copy.deepcopy(e)
                path_index = len(self.paths)
                self.paths.append(path)
                contained = np.where(self.overlap[:,e.index]==2)[0]
                for i in contained:
                    path.includes.add(i)
                    part = self.elements[i]
                    part.assigned_to.append(path_index)
                    self.assignments[i] += 1
        
        if len(self.SP)+len(self.EP)+len(self.SM)+len(self.EM) == 0:
            self.no_ends = True
        else:
            self.no_ends = False
        
        if len(self.paths) > 0:
            self.assign_weights()
    
    cpdef void assign_weights(self):
        """One round of Expectation Maximization: 
        Given existing weights of Paths and all Element assignments,
        Set Path weights as sum of assigned Elements * proportion of
        all weights of Paths that Element is assigned to.
        Runs once at initialization and once after each round of find_optimal_path()."""
        cdef int number_of_sources
        cdef Py_ssize_t i, p
        cdef np.ndarray priors, path_covs, sample_totals, proportions, cov_proportions
        cdef Element path, element
        number_of_sources = self.elements[0].weights.shape[0]
        priors = np.zeros(shape=(len(self.paths), number_of_sources))
        for i in range(len(self.paths)):
            priors[i,:] = self.paths[i].weights
        
        path_covs = np.sum(priors, axis=1, keepdims=True)
        cov_proportions = path_covs/np.sum(path_covs)
        sample_totals = np.sum(priors, axis=0)
        proportions = np.full(shape=(len(self.paths), number_of_sources), fill_value=cov_proportions)
        for i in np.where(sample_totals > 0)[0]:
            proportions[:,i] = priors[:,i]/sample_totals[i]
        
        for p in range(len(self.paths)):
            path = self.paths[p]
            path.weights = np.zeros(shape=(number_of_sources), dtype=np.float32)
            for i in path.includes:
                element = self.elements[i]
                path.weights += element.weights * proportions[p,:] * element.length
            
            path.weights /= path.length
    
    cpdef void assemble(self, float minimum_proportion):
        """Iteratively perform find_optimal_path() on the graph
        until the number of novel reads fails to exceed minimum_proportion
        of the reads at the locus. If minimum_proportion == 0, assemble()
        only terminates when every read is in a path."""
        cdef float threshold, total_bases_assigned
        cdef Element path
        
        total_bases_assigned = sum([path.bases for path in self.paths])
        path = self.find_optimal_path()
        if path is self.emptyPath:
            return
        
        threshold = self.bases*(1-minimum_proportion)
        total_bases_assigned += self.add_path(path)
        while total_bases_assigned < threshold:
            path = self.find_optimal_path()
            if path is self.emptyPath:
                total_bases_assigned = threshold
            else:
                total_bases_assigned += self.add_path(path)

    cpdef np.ndarray available_proportion(self, np.ndarray weights, Element element):
        """Given a path that wants to merge with the indexed element,
        calculate how much coverage is actually available to the path."""
        # Get the total cov of all already assigned paths
        cdef:
            np.ndarray assigned_weights, proportion
            int i
        if len(element.assigned_to) == 0: # No competition, all reads are available
            return element.all
        
        assigned_weights = np.copy(weights)
        proportion = np.ones(weights.shape[0], dtype=np.float32)
        for i in element.assigned_to:
            assigned_weights += self.paths[i].weights
        
        for i in np.where(assigned_weights > weights)[0]:
            proportion[i] = weights[i]/assigned_weights[i]
        
        return proportion
    
    cpdef void extend_path(self, Element path, tuple extension):
        """Merges the proper 
        """
        cdef Element extpath
        cdef int i
        cdef np.ndarray prior_weights, proportion
        prior_weights = np.copy(path.weights)
        for i in range(len(extension)):
            extpath = self.elements[extension[i]]
            proportion = self.available_proportion(prior_weights, extpath)
            path.merge(extpath, proportion)
    
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
    
    cpdef float dead_end(self, Element path, tuple extension):
        """Returns a multiplier that indicates how many termini can be
        reached by extending the path through extension:
        1 = neither end can be reached
        10 = one end can be reached
        100 = both ends can be reached."""
        cdef Element element
        cdef int i, strand
        cdef bint s_tag, e_tag
        cdef set includes, excludes, starts, ends
        if self.no_ends:
            return 1
        
        s_tag = path.s_tag
        e_tag = path.e_tag
        strand = path.strand
        excludes = copy.copy(path.excludes)
        for i in extension:
            element = self.elements[i]
            excludes.update(element.excludes)
            s_tag = s_tag or element.s_tag
            e_tag = e_tag or element.e_tag
            if strand == 0:
                strand = element.strand
        
        if s_tag and e_tag: # Both ends are already found
            return 100
        
        if not s_tag: # Try to extend from one of the known s_tag elements to path
            starts = [self.SM, self.SM|self.SP, self.SP][strand].difference(excludes)
            if len(starts) > 0: # At least one start is not excluded from the path
                s_tag = True
        
        if not e_tag:
            ends = [self.EM, self.EM|self.EP, self.EP][strand].difference(excludes)
            if len(ends) > 0: # At least one end is not excluded from the path
                e_tag = True
        
        return 1 * [1,10][s_tag] * [1,10][e_tag]
    
    cpdef list generate_extensions(self, Element path):
        """Defines all combinations of mutally compatible elements in
        the path's ingroup/outgroup that should be evaluated. Requires
        at least one each from ingroup and outgroup if they are nonempty."""
        cdef:
            np.ndarray ingroup, outgroup
            list pairs
            set ext_accounts, ext_members, ext_nonmembers, exclude, extensions, contained
            int i, o, c
            Element e, e_in, e_out, e_con
            (int, int) pair
        ingroup = np.array(sorted(path.ingroup))
        outgroup = np.array(sorted(path.outgroup))
        if len(ingroup) > 0:
            if len(outgroup) > 0:
                pairs = list(set([(i,o) for o in outgroup for i in ingroup if self.overlap[i,o] > -1]))
            else: # No outgroups, use path.index as other end of pair
                pairs = [(i,path.index) for i in ingroup]
        else: # No ingroups, use path.index as other end of pair
            pairs = [(path.index,o) for o in outgroup]
        
        # Make an extension set out of each pair by adding all elements contained by path+pair
        extensions = set()
        for pair in pairs:
            e_in = self.elements[pair[0]]
            e_out = self.elements[pair[1]]
            contained = e_in.outgroup | e_out.ingroup # Potential set of elements contained in the extension
            # Filter 1: All elements already included or excluded in the extension itself
            ext_accounts = e_in.includes | path.includes | e_out.includes | e_in.excludes | path.excludes | e_out.excludes
            contained.difference_update(ext_accounts)
            # Filter 2: All elements in the set that add information not contained in the extension
            stranded = e_in.strand != 0 or e_out.strand !=0 or path.strand != 0
            ext_members = e_in.members | path.members | e_out.members
            ext_nonmembers = e_in.nonmembers | path.nonmembers | e_out.nonmembers
            exclude = set([path.index])
            for c in contained:
                e_con = self.elements[c]
                if not stranded and e_con.strand != 0:
                    exclude.add(c)
                
                if not e_con.members.issubset(ext_members) or not e_con.nonmembers.issubset(ext_nonmembers):
                    exclude.add(c)
            
            contained.update([pair[0], pair[1]])
            contained.difference_update(exclude)
            extensions.add(tuple(sorted(list(contained))))
        
        return sorted(list(extensions))
    
    cpdef tuple best_extension(self, Element path, list extensions):
        cdef tuple ext, best_ext
        cdef float score, best_score
        best_ext = ()
        best_score = -1
        for ext in extensions:
            score = self.calculate_extension_score(path, ext)
            if score > best_score:
                best_ext = ext
                best_score = score
        
        return best_ext
    
    cpdef float calculate_extension_score(self, Element path, tuple extension):
        """Given a path and a set of Elements to extend from it, calculate the
        new weights of the extended path and return a score 
        """
        cdef:
            Element element
            int i
            set new_members
            float bases, new_cov, similarity, e_bases
            np.ndarray e_prop, e_weights, proportions, path_proportions
        new_members = set()
        bases = path.bases
        path_proportions = path.weights/path.cov
        proportions = path_proportions*path.bases
        for i in extension:
            element = self.elements[i]
            e_prop = self.available_proportion(path.weights, element)
            e_weights = e_prop*element.weights
            e_bases = np.sum(e_weights)*element.length
            bases += e_bases
            proportions += e_weights*element.length
            new_members.update(element.members.difference(path.members))
        
        proportions /= bases
        new_length = sum([path.frag_len[i] for i in new_members])
        # Calculate the new coverage (reads/base) of the extended path
        new_cov = bases/(path.length+new_length)
        similarity = 2 - np.sum(np.abs(path_proportions - proportions))
        dead_end_penalty = self.dead_end(path, extension)
        return new_cov * similarity * dead_end_penalty
    
    cpdef Element find_optimal_path(self, bint verbose=False):
        """Traverses the path in a greedy fashion from the heaviest element."""
        cdef Element currentPath, e
        cdef tuple ext
        cdef list extensions
        # Get the current working path (heaviest unassigned Element)
        currentPath = self.get_heaviest_element()
        if currentPath is self.emptyPath:
            return currentPath
        
        extensions = self.generate_extensions(currentPath)
        while len(extensions) > 0: # Extend as long as possible
            if len(extensions) == 1: # Only one option, do not evaluate
                self.extend_path(currentPath, extensions[0])
            else:
                ext = self.best_extension(currentPath, extensions)
                if verbose:
                    print("{} + {}".format(currentPath, ext))
                
                self.extend_path(currentPath, ext)
            
            extensions = self.generate_extensions(currentPath)
        
        return currentPath
    
    cpdef float add_path(self, Element path):
        """Evaluate what proportion of the compatible reads should be """
        cdef int i
        cdef float novel_bases = 0
        # Assign each included element to the path
        for i in path.includes:
            if self.assignments[i] == 0:
                novel_bases += self.elements[i].bases
            
            self.assignments[i] += 1
            self.elements[i].assigned_to.append(len(self.paths))
            self.elements[i].update()
        
        # Add the new path to the list of paths
        self.paths.append(path)
        self.assign_weights()
        return novel_bases
    
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
        
        self.assign_weights()


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

cdef class Element:
    """Represents a read or collection of reads in a Locus."""
    cdef public int index, length, IC, maxIC, left, right, number_of_elements, LM, RM
    cdef public list assigned_to
    cdef public char strand
    cdef public float cov, coverage, bases
    cdef public set members, nonmembers, ingroup, outgroup, excludes, includes, end_indices
    cdef public np.ndarray frag_len, weights, all
    cdef public bint complete, s_tag, e_tag, empty, is_spliced, has_gaps
    def __init__(self, int index, np.ndarray weights, char strand, np.ndarray membership, np.ndarray overlap, np.ndarray frag_len, int maxIC):
        cdef Py_ssize_t i
        cdef char m, overOut, overIn
        self.is_spliced = False                       # Default: the path has no discontinuities
        self.index = self.left = self.right = index   # Initialize left, right, and index
        self.number_of_elements = overlap.shape[0]    # Total number of nodes in the graph
        self.frag_len = frag_len                      # Length of the fragment is provided
        self.includes = set([self.index])             # Which Elements are part of this Element
        self.excludes = set()                         # Which Elements are incompatible with this Element
        self.weights = weights                        # Array of read coverage per Source
        self.strand = strand                          # +1, -1, or 0 to indicate strand of path
        self.length = 0                               # Number of nucleotides in the path
        self.assigned_to = []                         # List of Path indices this Element is a part of
        self.complete = False                         # Represents an entire end-to-end transcript
        self.has_gaps = False                         # Is missing information
        self.members = set()                          # Set of Member indices contained in this Element
        self.nonmembers = set()                       # Set of Members indices incompatible with this Element
        self.ingroup = set()                          # Set of compatible upstream Elements
        self.outgroup = set()                         # Set of Compatible downstream Elements
        self.all = np.ones(shape=self.weights.shape[0], dtype=np.float32)
        if index == -1:                               # Special Element emptyPath: placeholder for null values
            self.empty = True
            self.maxIC = 0
            self.end_indices = set()
        else:
            self.empty = False
            self.maxIC = maxIC
            self.end_indices = set(range(self.maxIC-4, self.maxIC))
            for i in range(self.maxIC):
                m = membership[i]
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
    
    def __eq__(self, other): return self.bases == other.bases
    def __ne__(self, other): return self.bases != other.bases
    def __gt__(self, other): return self.bases >  other.bases
    def __ge__(self, other): return self.bases >= other.bases
    def __lt__(self, other): return self.bases <  other.bases
    def __le__(self, other): return self.bases <= other.bases

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
        
        summed_element.merge(other, self.all)
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
        
        self.cov = sum(self.weights)
        self.bases = self.cov*self.length
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
    
    cpdef set uniqueNonmembers(self, Element other):
        """Given a second Element, return a set of frags that
        are only in self and not in other"""
        return self.nonmembers.difference(other.nonmembers)
    
    cpdef int uniqueInformation(self, Element other):
        """Given a second Element, return a set of frags that
        are only in self and not in other"""
        return len(self.uniqueMembers(other)|self.uniqueNonmembers(other))
    
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
    
    cpdef void merge(self, Element other, np.ndarray proportion):
        """Add an Element to this one, combining their membership and reads
        through in-place updates of self."""
        cdef set covered, unique
        cdef bint extendedLeft, extendedRight
        cdef int o, i, leftMember, rightMember
        cdef float old_length
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
        
        old_length = self.length
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
        
        self.weights = (other.weights*other.length*proportion + self.weights*old_length)/self.length
        self.update()
        if self.strand == 1: # Enforce directionality of edges
            self.outgroup = set([o for o in self.outgroup if o > self.right])
            self.ingroup = set([i for i in self.ingroup if i < self.left])
        elif self.strand == -1:
            self.outgroup = set([o for o in self.outgroup if o < self.left])
            self.ingroup = set([i for i in self.ingroup if i > self.right])
