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
    cdef public np.ndarray assignments, overlap, end_reachability
    cdef readonly int number_of_elements, maxIC
    cdef Element emptyPath
    cdef public float bases, input_bases, dead_end_penalty
    cdef public set SP, SM, EP, EM
    cdef public bint no_ends, naive, partial_coverage
    def __init__(self, np.ndarray overlap_matrix, np.ndarray membership_matrix, source_weight_array, member_weight_array, strands, lengths, naive=False, dead_end_penalty=0.1, partial_coverage=True):
        """Constructs a forward and reverse directed graph from the
        connection values (ones) in the overlap matrix.
        Additionally, stores the set of excluded edges for each node as an 'antigraph'
        """
        cdef int e_index, path_index, i
        cdef Element e, path, part
        self.emptyPath = Element(-1, np.array([0]), np.array([0]), 0, np.array([0]), np.array([0]), np.array([0]), 0)
        self.dead_end_penalty = dead_end_penalty
        self.SP, self.SM, self.EP, self.EM = set(), set(), set(), set()
        self.overlap = overlap_matrix
        self.number_of_elements = self.overlap.shape[0]
        self.maxIC = membership_matrix.shape[1]
        self.naive = naive
        self.partial_coverage = partial_coverage
        if self.naive:
                source_weight_array = np.sum(source_weight_array, axis=1, keepdims=True)
        
        self.elements = [Element(
            i, source_weight_array[i,:], member_weight_array[i,:], strands[i],
            membership_matrix[i,:], self.overlap, lengths, self.maxIC
        ) for i in range(self.number_of_elements)] # Generate an array of Element objects
        self.assignments = np.zeros(shape=self.number_of_elements, dtype=np.int32)
        self.paths = []
        self.input_bases = sum([e.bases for e in self.elements])
        self.check_for_full_paths()
        self.penalize_dead_ends()
        self.bases = sum([e.bases for e in self.elements])
        if self.bases > 0:
            self.resolve_containment()
    
    cpdef void penalize_dead_ends(self):
        """Perform a breadth-first search from all starts and all ends.
        The weight of all elements unreachable by each search is multiplied
        by the dead_end_penalty, for a maximum penalty of dead_end_penalty^2"""
        cdef Element element
        cdef np.ndarray reached_from_start, reached_from_end
        cdef float original_bases
        cdef int i, Sp, Ep, Sm, Em
        cdef set strand
        self.end_reachability = np.zeros(shape=(4,len(self.elements)), dtype=np.bool)
        Sp, Ep, Sm, Em = 0, 0, 0, 0
        queue = deque(maxlen=len(self.elements))
        strand = set([0,1])
        queue.extend(sorted(self.SP))
        while queue:
            v = queue.popleft()
            self.end_reachability[Sp, v] = True
            self.end_reachability[Sp, sorted(self.elements[v].contains)] = True
            for w in self.elements[v].outgroup:
                if not self.end_reachability[Sp, w]:
                    if self.elements[w].strand in strand:
                        self.end_reachability[Sp, w] = True
                        queue.append(w)
        
        queue.clear()
        queue.extend(self.ends_plus)
        while queue:
            v = queue.popleft()
            self.end_reachability[Ep, v] = True
            self.end_reachability[Ep, sorted(self.elements[v].contains)] = True
            for w in self.elements[v].ingroup:
                if not self.end_reachability[Ep, w]:
                    if self.elements[w].strand in strand:
                        self.end_reachability[Ep, w] = True
                        queue.append(w)
        
        queue.clear()
        queue.extend(self.starts_minus)
        strand = set([0,-1])
        while queue:
            v = queue.popleft()
            self.end_reachability[Sm, v] = True
            self.end_reachability[Sm, sorted(self.elements[v].contains)] = True
            for w in self.elements[v].ingroup:
                if not self.end_reachability[Sm, w]:
                    if self.elements[w].strand in strand:
                        self.end_reachability[Sm, w] = True
                        queue.append(w)
        
        queue.clear()
        queue.extend(self.ends_minus)
        while queue:
            v = queue.popleft()
            self.end_reachability[Em, v] = True
            self.end_reachability[Em, sorted(self.elements[v].contains)] = True
            for w in self.elements[v].outgroup:
                if not self.end_reachability[Em, w]:
                    if self.elements[w].strand in strand:
                        self.end_reachability[Em, w] = True
                        queue.append(w)
        
        reached_from_start = np.logical_or(self.end_reachability[Sp,:] , self.end_reachability[Sm,:])
        reached_from_end = np.logical_or(self.end_reachability[Ep,:], self.end_reachability[Em,:])
        for i in range(len(self.elements)):
            element = self.elements[i]
            if not reached_from_start[i]:
                if self.dead_end_penalty == 0:
                    self.zero_element(i)
                else:
                    element.source_weights *= self.dead_end_penalty
                    element.member_weights *= self.dead_end_penalty
                    element.cov *= self.dead_end_penalty
                    original_bases = element.bases
                    element.bases *= self.dead_end_penalty
                    self.bases -= original_bases-element.bases
            
            if not reached_from_end[i]:
                if self.dead_end_penalty == 0:
                    self.zero_element(i)
                else:
                    element.source_weights *= self.dead_end_penalty
                    element.member_weights *= self.dead_end_penalty
                    element.cov *= self.dead_end_penalty
                    original_bases = element.bases
                    element.bases *= self.dead_end_penalty
                    self.bases -= original_bases-element.bases
    
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
    
    cpdef void resolve_containment(self):
        """Given a overlap matrix, 'bubble up' the weight of
        all reads that have one or more 'contained by' relationships
        to other reads. Pass from highest complexity reads down, assigning
        weight proportional to the existing weight.
        The resulting matrix should contain only overlaps, exclusions, and unknowns."""
        cdef:
            np.ndarray default_proportions, proportions
            list contained, resolve_order, container_indices, maxmember
            set zeros, containers, incompatible
            float total_container_weight, m
            Element element
            Py_ssize_t resolve, i, c
            (int, int, int) sorttuple
        
        zeros = set([element.index for element in self.elements if element.cov == 0])
        contained = [i for i in range(self.number_of_elements) if len(self.elements[i].contained)>0] # Identify reads that are contained
        resolve_order = [sorttuple[2] for sorttuple in sorted([(-self.elements[c].IC, len(self.elements[c].contained), c) for c in contained])] # Rank them by decreasing number of members
        for resolve in resolve_order:
            element = self.elements[resolve]
            if element.cov == 0:continue
            containers = element.contained.difference(zeros)
            # Get the set of reads incompatible with all containers but that do not exclude i
            incompatible =  set(range(self.number_of_elements))
            for c in containers:
                incompatible.intersection_update(self.elements[c].excludes)
            
            incompatible.difference_update(element.excludes|zeros)
            if len(incompatible) == 0: # Special case, all weight goes to containers
                if len(containers) == 1:
                    self.elements[containers.pop()].merge(element, element.all)
                else: # Evaluate how much weight goes to each container
                    container_indices = sorted(containers)
                    maxmember = [self.elements[c].cov for c in container_indices]
                    total_container_weight = sum(maxmember)
                    default_proportions = np.array([m/total_container_weight for m in maxmember], dtype=np.float32)
                    proportions = np.zeros(shape=(len(containers), element.all.shape[0]), dtype=np.float32)
                    for i in range(len(container_indices)):
                        c = container_indices[i]
                        proportions[i,:] = self.normalize(self.elements[c].source_weights)*default_proportions[i]
                    
                    for i in range(element.all.shape[0]):
                        if np.sum(proportions[:,i]) == 0:
                            proportions[:,i] = default_proportions
                        else:
                            proportions[:,i] = self.normalize(proportions[:,i])
                    
                    for i in range(len(container_indices)):
                        self.elements[container_indices[i]].merge(element, proportions[i,:])
                
                self.zero_element(resolve)
                zeros.add(resolve)
    
    cpdef void zero_element(self, int index):
        """Given an element's index, remove all references to it without
        deleting it from the list of elements."""
        cdef Element element
        self.assignments[index] = -1
        element = self.elements[index]
        element.cov = 0
        element.bases = 0
        element.source_weights -= element.source_weights
        element.member_weights -= element.member_weights
        for e in self.elements:
            e.ingroup.discard(index)
            e.outgroup.discard(index)
            e.contains.discard(index)
            e.contained.discard(index)
            e.includes.discard(index)
            e.excludes.discard(index)
    
    cpdef void assign_weights(self):
        """One round of Expectation Maximization: 
        Given existing weights of Paths and all Element assignments,
        Set Path weights as sum of assigned Elements * proportion of
        all weights of Paths that Element is assigned to.
        Runs once at initialization and once after each round of find_optimal_path()."""
        cdef int number_of_sources
        cdef Py_ssize_t i, p, j
        cdef np.ndarray priors, path_covs, sample_totals, proportions, cov_proportions, assignment_proportions
        cdef Element path, element
        number_of_sources = self.elements[0].source_weights.shape[0]
        priors = np.zeros(shape=(len(self.paths), number_of_sources))
        for i in range(len(self.paths)):
            priors[i,:] = self.paths[i].source_weights
            self.paths[i].source_weights = np.zeros(shape=(number_of_sources), dtype=np.float32)
        
        path_covs = np.sum(priors, axis=1, keepdims=True)
        cov_proportions = path_covs/np.sum(path_covs)
        sample_totals = np.sum(priors, axis=0)
        proportions = np.full(shape=(len(self.paths), number_of_sources), fill_value=cov_proportions)
        for i in np.where(sample_totals > 0)[0]:
            proportions[:,i] = priors[:,i]/sample_totals[i]
        
        for i in np.argsort(self.assignments): # Assign FULL weight of each assigned element
            if self.assignments[i] == 0: continue
            element = self.elements[i]
            if self.assignments[i] == 1: # No proportions needed, assign all weight to 1 path
                self.paths[element.assigned_to[0]].source_weights += element.source_weights * element.length
            else: # Assigned paths must compete for element's source_weights
                assignment_proportions = proportions[element.assigned_to,:]
                assignment_proportions = np.apply_along_axis(self.normalize, 0, assignment_proportions)
                for j in range(len(element.assigned_to)):
                    path = self.paths[element.assigned_to[j]]
                    path.source_weights += element.source_weights * element.length * assignment_proportions[j,:]
        
        for path in self.paths: # Update path source_weights
            path.source_weights /= path.length
            path.bases = sum(path.source_weights)*path.length
    
    cpdef void assemble(self, float minimum_proportion):
        """Iteratively perform find_optimal_path() on the graph
        until the number of novel reads fails to exceed minimum_proportion
        of the reads at the locus. If minimum_proportion == 0, assemble()
        only terminates when every read is in a path."""
        cdef float threshold, total_bases_assigned, novel_bases
        cdef Element path
        
        total_bases_assigned = sum([self.elements[i].bases for i in np.where(self.assignments>0)[0]])
        threshold = self.bases*(1-minimum_proportion)
        while total_bases_assigned < threshold:
            path = self.find_optimal_path(minimum_proportion)
            if path is self.emptyPath:
                total_bases_assigned = threshold
            else:
                novel_bases = self.add_path(path)
                if novel_bases == 0:
                    total_bases_assigned = threshold
                else:
                    total_bases_assigned += novel_bases
    
    cpdef np.ndarray available_proportion(self, np.ndarray weights, Element element):
        """Given a path that wants to merge with the indexed element,
        calculate how much coverage is actually available to the path."""
        # Get the total cov of all already assigned paths
        cdef:
            np.ndarray assigned_weights, proportion
            int i
            float free_weight, min_weight, coverage_over_element, coverage_outside_element
            Element path
            # set outside_element
        if len(element.assigned_to) == 0: # No competition, all reads are available
            return element.all
        
        assigned_weights = np.copy(weights)
        for i in element.assigned_to:
            path = self.paths[i]
            assigned_weights += path.source_weights
            # if len(element.members.intersection(path.bottleneck)) > 0: # Element is in path's bottleneck
            #     return np.zeros(element.all.shape[0], dtype=np.float32)
        
        proportion = np.ones(weights.shape[0], dtype=np.float32)
        for i in np.where(assigned_weights > weights)[0]:
            proportion[i] = weights[i]/assigned_weights[i]
        
        return proportion
    
    cpdef float available_bases(self, np.ndarray weights, Element element):
        """Given a path to merge, calculate the number of bases available for merging"""
        cdef np.ndarray proportion
        proportion = self.available_proportion(weights, element)
        return np.sum(element.source_weights*proportion)*element.length
    
    cpdef void extend_path(self, Element path, tuple extension):
        """Merges the proper 
        """
        cdef Element extpath
        cdef int i
        cdef np.ndarray prior_weights, proportion
        prior_weights = np.copy(path.source_weights)
        for i in range(len(extension)):
            extpath = self.elements[extension[i]]
            proportion = self.available_proportion(prior_weights, extpath)
            path.merge(extpath, proportion)
    
    cpdef Element get_heaviest_element(self):
        cdef Element best_element, new_element
        cdef float most_bases
        cdef int c
        cdef np.ndarray available_elements = np.where(self.assignments==0)[0]
        if len(available_elements) == 0:
            return self.emptyPath
        
        best_element = self.elements[available_elements[0]]
        most_cov = best_element.bases
        for c in best_element.contains:
            if c != best_element.index:
                most_cov += self.available_bases(best_element.source_weights, self.elements[c])
        
        most_cov /= best_element.length
        for i in available_elements:
            new_element = self.elements[i]
            new_cov = new_element.bases
            for c in new_element.contains:
                if c != new_element.index:
                    new_cov += self.available_bases(new_element.source_weights, self.elements[c])
            
            new_cov /= new_element.length
            # print(i, new_element, new_cov)
            if new_cov > most_cov:
                best_element = new_element
                most_cov = new_cov
            elif new_cov == most_cov: # Break ties by complexity
                if new_element.IC > best_element.IC:
                    best_element = new_element
                    most_cov = new_cov
        
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
        for i in extension:
            element = self.elements[i]
            s_tag = s_tag or element.s_tag
            e_tag = e_tag or element.e_tag
            if strand == 0:
                strand = element.strand
        
        if s_tag and e_tag: # Both ends are already found
            return 1
        
        if not s_tag: # Check that BFS from the + and/or - Start reached all extension indices
            if strand >= 0:
                s_tag = np.logical_and(self.end_reachability[0,extension])
            
            if strand <= 0 and not s_tag:
                s_tag = np.logical_and(self.end_reachability[2,extension])
        
        if not e_tag:
            if strand >= 0:
                e_tag = np.logical_and(self.end_reachability[1,extension]) 
            
            if strand <= 0 and not e_tag:
                e_tag = np.logical_and(self.end_reachability[3,extension])
        
        return [self.dead_end_penalty,1.][s_tag] * [self.dead_end_penalty,1.][e_tag]
    
    cpdef list generate_extensions(self, Element path):
        """Defines all combinations of mutally compatible elements in
        the path's ingroup/outgroup that should be evaluated. Requires
        at least one each from ingroup and outgroup if they are nonempty."""
        cdef:
            np.ndarray ingroup, outgroup
            list pairs
            set ext_accounts, ext_members, ext_nonmembers, exclude, contained
            int i, o, c
            Element e, e_in, e_out, e_con
            (int, int) pair
            tuple freebies, ext
            dict extdict
            str exthash
        extdict = {}
        ingroup = np.array(sorted(path.ingroup))
        outgroup = np.array(sorted(path.outgroup))
        freebies = tuple(path.contains.difference(path.includes))
        if len(freebies) > 0:
            self.extend_path(path, freebies)
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
        for pair in pairs:
            e_in = self.elements[pair[0]]
            e_out = self.elements[pair[1]]
            contained = e_in.outgroup | e_out.ingroup | e_in.contains | e_out.contains # Potential set of elements contained in the extension
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
            ext = tuple(sorted(list(contained)))
            if len(ext) == 0 or ext_members.issubset(path.members):continue
            if len(ext) > 1 or len(self.elements[ext[0]].uniqueMembers(path)) > 0:
                exthash = '_'.join([','.join([str(i) for i in sorted(ext_members)]), ','.join([str(i) for i in sorted(ext_nonmembers)])])
                if len(ext) > len(extdict.get(exthash, ())):
                    extdict[exthash] = ext
        
        # Final check: If >0 extensions go both ways, remove the extensions that don't
        extensions = sorted([(sum([self.elements[i].cov for i in ext]),ext) for ext in list(extdict.values())],reverse=True)
        return extensions
    
    cpdef tuple best_extension(self, Element path, list extensions, float minimum_proportion):
        cdef tuple ext, best_ext
        cdef float score, best_score
        best_ext = ()
        best_score = 0
        for cov,ext in extensions:
            if cov >= best_score:
                score = self.calculate_extension_score(path, ext, minimum_proportion)
                if score > best_score or (score == best_score and len(ext) > len(best_ext)):
                    best_ext = ext
                    best_score = score
            else:
                break
            
        return best_ext
    
    # cpdef float calculate_extension_score(self, Element path, tuple extension, float minimum_proportion, float prune):
    #     """Given a path and a set of Elements to extend from it, calculate the
    #     new weights of the extended path and return a score 
    #     """
    #     cdef:
    #         Element element
    #         int i, outgroup_bases, excluded
    #         set new_members, extension_outgroup, excluded_members, extension_excludes, new_covered_indices
    #         float bases, new_bases, extension_bases, excluded_bases, score, source_similarity, e_cov, e_bases, ext_cov, ext_jcov, path_jcov, junction_delta, dead_end_penalty, excluded_cov
    #         np.ndarray e_prop, e_weights, proportions, path_proportions, new_member_weights, combined_member_coverage
    #         list shared_members
    #     new_members = set()
    #     bases = path.bases
    #     path_proportions = path.source_weights/path.cov
    #     proportions = path_proportions*path.bases
    #     new_member_weights = np.zeros(path.member_weights.shape[0], dtype=np.float32)
    #     extension_bases = 0
    #     extension_outgroup = set()
    #     extension_excludes = set()
    #     new_covered_indices = set()
    #     for i in extension:
    #         element = self.elements[i]
    #         new_covered_indices.update(element.covered_indices)
    #         extension_outgroup.update((element.outgroup|element.ingroup).difference(path.excludes|path.includes))
    #         extension_excludes.update(element.excludes)
    #         e_prop = self.available_proportion(path.source_weights, element)
    #         e_weights = e_prop*element.source_weights
    #         e_cov = np.sum(e_weights)
    #         e_bases = e_cov*element.length
    #         extension_bases += element.bases
    #         bases += e_bases
    #         new_member_weights += element.member_weights*e_cov/element.cov
    #         proportions += e_weights*element.length
    #         new_members.update(element.members.difference(path.members))
        
    #     proportions /= bases
    #     new_length = sum([path.frag_len[i] for i in new_members])
    #     # Calculate the new coverage (reads/base) of the extended path
    #     new_bases = bases - path.bases
    #     if self.partial_coverage: # Include overlapping weight of all members compatible with path + extension
    #         extension_outgroup.difference_update(set(extension)|extension_excludes)
    #         for i in extension_outgroup: # Add bases of the overlapping portions of all compatible outgroups
    #             element = self.elements[i]
    #             e_prop = self.available_proportion(path.source_weights, element)
    #             shared_members = sorted(element.members.intersection(new_members|path.members))
    #             shared_length = np.sum(path.frag_len[shared_members])
    #             outgroup_bases = np.sum(element.source_weights*e_prop)*shared_length
    #             new_bases += outgroup_bases
    #             extension_bases += np.sum(element.source_weights)*shared_length
    #             new_member_weights[shared_members] += element.member_weights[shared_members]*np.sum(element.source_weights*e_prop)
        
    #     if new_bases < minimum_proportion * extension_bases:
    #         return 0
        
    #     ext_cov = new_bases / new_length
    #     extension_excludes.difference_update(path.excludes)
    #     if len(extension_excludes) > 0:
    #         excluded_bases = 0
    #         excluded_members = set()
    #         for excluded in extension_excludes:
    #             element = self.elements[excluded]
    #             excluded_bases += self.available_bases(path.source_weights, element)
    #             excluded_members.update(element.members)
        
    #         excluded_length = np.sum(path.frag_len[sorted(excluded_members)])
    #         excluded_cov = excluded_bases/excluded_length
    #         exclusion_penalty = ext_cov/(ext_cov+excluded_cov)
    #     else:
    #         excluded_cov = 0
    #         exclusion_penalty = 1
        
    #     combined_member_coverage = np.add(path.member_weights,new_member_weights)[sorted(path.covered_indices.union(new_covered_indices))]
    #     weakest_link_penalty = np.min(combined_member_coverage)/np.mean(combined_member_coverage)
    #     source_similarity = 2 - np.sum(np.abs(path_proportions - proportions))
    #     dead_end_penalty = self.dead_end(path, extension)
    #     score = ext_cov * source_similarity * weakest_link_penalty * dead_end_penalty * exclusion_penalty
    #     return score
    
    cpdef float calculate_extension_score(self, Element path, tuple extension, float minimum_proportion):
        """Given a path and a set of Elements to extend from it, calculate the
        new weights of the extended path and return a score 
        """
        cdef:
            Element element
            int i
            set extension_excludes, new_covered_indices
            float div, score, source_similarity, ext_cov, dead_end_penalty, variance_penalty
            np.ndarray ext_proportions, e_prop, path_proportions, combined_member_coverage
            list shared_members, excluded_cov
        if len(extension)==0:return 0
        ext_member_weights = np.zeros(path.member_weights.shape[0], dtype=np.float32)
        ext_proportions = np.zeros(path.source_weights.shape[0], dtype=np.float32)
        new_covered_indices = set()
        extension_excludes = set()
        div = 1/len(extension)
        for i in extension:
            element = self.elements[i]
            new_covered_indices.update(element.covered_indices)
            extension_excludes.update(element.excludes)
            e_prop = self.available_proportion(path.source_weights, element)
            ext_proportions += self.normalize(e_prop*element.source_weights)*div
            ext_member_weights += element.member_weights*np.sum(ext_proportions)
        
        ext_cov = np.max(ext_member_weights[sorted(new_covered_indices)])
        extension_excludes.difference_update(path.excludes)
        if len(extension_excludes) > 0:
            excluded_cov = [self.elements[i].cov for i in extension_excludes]
            exclusion_penalty = ext_cov/(ext_cov+sum(excluded_cov))
        else:
            exclusion_penalty = 1
        
        combined_member_coverage = np.add(path.member_weights,ext_member_weights)[sorted(path.covered_indices.union(new_covered_indices))]
        variance_penalty = np.mean(combined_member_coverage)/np.max(combined_member_coverage)
        path_proportions = self.normalize(path.source_weights)
        source_similarity = .5*(2 - np.sum(np.abs(path_proportions - ext_proportions)))
        dead_end_penalty = self.dead_end(path, extension)
        score = ext_cov * source_similarity * weakest_link_penalty * dead_end_penalty * exclusion_penalty
        return score
    
    cpdef Element find_optimal_path(self, float minimum_proportion, bint verbose=False):
        """Traverses the path in a greedy fashion from the heaviest element."""
        cdef Element currentPath, e
        cdef tuple ext
        cdef list extensions
        # Get the current working path (heaviest unassigned Element)
        currentPath = self.get_heaviest_element()
        self.extend_path(currentPath, tuple(sorted(currentPath.contains)))
        
        if currentPath is self.emptyPath:
            return currentPath
        
        extensions = self.generate_extensions(currentPath)
        while len(extensions) > 0: # Extend as long as possible
            if len(extensions) == 1: # Only one option, do not evaluate
                if self.calculate_extension_score(currentPath, extensions[0][1], minimum_proportion) == 0:break
                self.extend_path(currentPath, extensions[0][1])
            else:
                ext = self.best_extension(currentPath, extensions, minimum_proportion)
                if verbose:print("{} + {}".format(currentPath, ext))
                if len(ext) == 0:break
                self.extend_path(currentPath, ext)
            
            extensions = self.generate_extensions(currentPath)
        
        if verbose:print(currentPath)
        # self.trim_ends(currentPath)
        return currentPath
    
    # cpdef void trim_ends(self, Element path):
    #     """If a path has malformed ends, check if it is possible to back up to a
    #     bypassed start/end site without crossing a splice junction."""
    #     pass
    #     cdef Element element, trim_element, start_element, end_element
    #     cdef np.ndarray members
    #     cdef int left_exon_border, right_exon_border, bypassed, m, lastm
    #     cdef set elements_to_trim
    #     cdef list candidates, repair_elements
    #     if path.complete or path.strand==0:return
    #     left_exon_border = -1
    #     right_exon_border = -1
    #     lastm = -1
    #     members = np.array(sorted(path.members.difference(path.end_indices)), dtype=np.int32)
    #     for m in members:
    #         if lastm == -1:lastm = m
    #         if m > lastm+1:
    #             if left_exon_border==-1:left_exon_border = lastm
    #             right_exon_border = m
            
    #         lastm = m
        
    #     if left_exon_border == -1:left_exon_border=path.RM # single-exon path
    #     if right_exon_border == -1:right_exon_border=path.LM # single-exon path
    #     start_to_repair = self.emptyPath
    #     if not path.s_tag: # Check if there is a bypassed start in the first exon
    #         if path.strand == 1:
    #             candidates = self.starts_plus
    #             for bypassed in candidates:
    #                 element = self.elements[bypassed]
    #                 if element.LM > left_exon_border:continue
    #                 if element.s_tag and element.LM > path.LM: # Start is in the bypassed range
    #                     if start_to_repair is self.emptyPath:
    #                         start_to_repair = element
    #                     elif element.LM < start_to_repair.LM: # New start is upstream of the start already found
    #                         start_to_repair = element
    #         elif path.strand == -1:
    #             candidates = self.starts_minus
    #             for bypassed in candidates:
    #                 element = self.elements[bypassed]
    #                 if element.RM < right_exon_border:continue
    #                 if element.s_tag and element.RM < path.RM:
    #                     if start_to_repair is self.emptyPath:
    #                         start_to_repair = element
    #                     elif element.RM > start_to_repair.RM: # New start is upstream of the start already found
    #                         start_to_repair = element
        
    #     end_to_repair = self.emptyPath
    #     if not path.e_tag: # Check if there is a bypassed end in the first exon
    #         if path.strand == -1:
    #             candidates = self.ends_minus
    #             for bypassed in candidates:
    #                 element = self.elements[bypassed]
    #                 if element.LM > left_exon_border:continue
    #                 if element.e_tag and element.LM > path.LM:
    #                     if end_to_repair is self.emptyPath:
    #                         end_to_repair = element
    #                     elif element.LM < end_to_repair.LM: # New end is downstream
    #                         end_to_repair = element
    #         elif path.strand == 1:
    #             candidates = self.ends_plus
    #             for bypassed in candidates:
    #                 element = self.elements[bypassed]
    #                 if element.RM < right_exon_border:continue
    #                 if element.e_tag and element.RM < path.RM:
    #                     if end_to_repair is self.emptyPath:
    #                         end_to_repair = element
    #                     elif element.RM > end_to_repair.RM: # New end is downstream
    #                         end_to_repair = element
        
    #     # Remove all included/excluded elements from path to make room for the repaired ends
    #     repair_elements = [element for element in [start_to_repair, end_to_repair] if element is not self.emptyPath]
    #     members_to_trim = set()
    #     nonmembers_to_trim = set()
    #     for element in repair_elements:
    #         elements_to_trim = path.includes.intersection(element.excludes)
    #         for e in elements_to_trim:
    #             trim_element = self.elements[e]
    #             members_to_trim.update(trim_element.members.intersection(element.nonmembers))
    #             nonmembers_to_trim.update(trim_element.nonmembers.intersection(element.members))
            
    #         path.members.difference_update(members_to_trim)
    #         path.nonmembers.difference_update(nonmembers_to_trim)
    #         path.includes.difference_update(elements_to_trim)
    #         path.contains.difference_update(elements_to_trim)
    #         path.excludes.remove(element.index)
    #         path.merge(element, self.available_proportion(path.source_weights, element))

    #     path.update()
    
    cpdef float add_path(self, Element path):
        """Evaluate what proportion of the compatible reads should be """
        cdef int i
        cdef float novel_bases = 0
        cdef Element existing_path
        # Assign each included element to the path
        for existing_path in self.paths:
            if path.compatible(existing_path):
                return 0
        
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
    
    cpdef np.ndarray normalize(self, np.ndarray arr):
        cdef float arrsum
        arrsum = np.sum(arr)
        if arrsum > 0:
            return arr/arrsum
        else:
            return arr

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
    cdef public float cov, bases, bottleneck_weight
    cdef public set members, nonmembers, ingroup, outgroup, contains, contained, excludes, includes, end_indices, covered_indices, bottleneck
    cdef public np.ndarray frag_len, source_weights, member_weights, all
    cdef public bint complete, s_tag, e_tag, empty, is_spliced, has_gaps
    def __init__(self, int index, np.ndarray source_weights, np.ndarray member_weights, char strand, np.ndarray membership, np.ndarray overlap, np.ndarray frag_len, int maxIC):
        cdef Py_ssize_t i
        cdef char m, overOut, overIn
        self.is_spliced = False                       # Default: the path has no discontinuities
        self.index = self.left = self.right = index   # Initialize left, right, and index
        self.number_of_elements = overlap.shape[0]    # Total number of nodes in the graph
        self.frag_len = frag_len                      # Length of the fragment is provided
        self.includes = set([self.index])             # Which Elements are part of this Element
        self.excludes = set()                         # Which Elements are incompatible with this Element
        self.source_weights = np.copy(source_weights) # Array of read coverage per Source
        self.member_weights = np.copy(member_weights) # Array of read coverage of all members
        self.strand = strand                          # +1, -1, or 0 to indicate strand of path
        self.length = 0                               # Number of nucleotides in the path
        self.assigned_to = []                         # List of Path indices this Element is a part of
        self.complete = False                         # Represents an entire end-to-end transcript
        self.has_gaps = False                         # Is missing information
        self.members = set()                          # Set of Member indices contained in this Element
        self.nonmembers = set()                       # Set of Members indices incompatible with this Element
        self.ingroup = set()                          # Set of compatible upstream Elements
        self.outgroup = set()                         # Set of Compatible downstream Elements
        self.contains = set()
        self.contained = set()
        self.all = np.ones(shape=self.source_weights.shape[0], dtype=np.float32)
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
            
            for i in range(self.number_of_elements):
                if i == self.index:
                    continue
                
                overOut = overlap[self.index, i]
                if overOut == -1:
                    self.excludes.add(i)
                    continue
                elif overOut >= 1:
                    self.outgroup.add(i)
                    if overOut == 2:
                        self.contained.add(i)
                
                overIn = overlap[i, self.index]
                if overIn >= 1:
                    self.ingroup.add(i)
                    if overIn == 2:
                        self.contains.add(i)
            
            self.update()
    
    def __repr__(self):
        chars = [' ']*self.maxIC
        strand = {-1:'-', 0:'.', 1:'+'}[self.strand]
        for m in self.members:
            chars[m] = '*'
        
        for n in self.nonmembers:
            chars[n] = '_'
        
        return '|{}| {}-{} ({})'.format(''.join(chars),self.left,self.right,strand)
    
    cdef str span_to_string(self, (int, int) span):
        """Converts a tuple of two ints to a string connected by ':'"""
        return '{}:{}'.format(span[0], span[1])
    
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
        
        summed_element.merge(other, self.all)
        return summed_element
    
    cpdef void update(self):
        cdef int n, lastn
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
        
        self.bases = sum(self.source_weights)*self.length
        self.IC = len(self.members) + len(self.nonmembers)
        if self.IC == self.maxIC:
            self.complete = True
        else:
            self.complete = False
        
        self.covered_indices = self.members.difference(self.end_indices) # Covered regions (excluding starts/ends)
        for n in sorted(self.nonmembers):
            lastn = -1
            if n > self.LM and n < self.RM:
                self.is_spliced = True
                if n > lastn+1: # Only add one representative nonmember per intron
                    self.covered_indices.add(n)
        
        self.cov = np.max(self.member_weights[sorted(self.covered_indices)])
        # self.bottleneck_weight = np.min(self.member_weights[sorted(self.covered_indices)])
        # self.bottleneck = set(np.where(self.member_weights == self.bottleneck_weight)[0])
    
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
        
        if self.strand == 0 and other.strand != 0:
            self.strand = other.strand
        
        old_length = self.length
        self.contains.update(other.contains)
        self.contained.intersection_update(other.contained)
        unique = other.uniqueMembers(self)
        # Update Membership
        self.nonmembers.update(other.nonmembers)
        # Update Overlaps
        self.excludes.update(other.excludes) # Sum of exclusions
        self.includes.update(other.includes) # Sum of inclusions
        self.contains.update(self.includes)
        if len(unique) > 0: 
            self.outgroup.update(other.outgroup.difference(self.excludes))
            self.ingroup.update(other.ingroup.difference(self.excludes))
            for f in unique: # Append each frag from other to self
                self.members.add(f)
                self.length += self.frag_len[f]
            
            # Update the left and right borders of the Element
            self.right = max(self.right, other.right)
            self.left = min(self.left, other.left)
        
        self.outgroup.difference_update(self.contains|self.excludes)
        self.ingroup.difference_update(self.contains|self.excludes)
        self.source_weights = (other.source_weights*other.length*proportion + self.source_weights*old_length)/self.length
        self.member_weights += other.member_weights*np.sum(other.source_weights*proportion)/np.sum(other.source_weights)
        self.right = max(self.right, other.right)
        self.left = min(self.left, other.left)
        self.update()
        # if self.strand == 1: # Enforce directionality of edges
        #     self.outgroup = set([o for o in self.outgroup if o > self.right])
        #     self.ingroup = set([i for i in self.ingroup if i < self.left])
        # elif self.strand == -1:
        #     self.outgroup = set([o for o in self.outgroup if o < self.left])
        #     self.ingroup = set([i for i in self.ingroup if i > self.right])


