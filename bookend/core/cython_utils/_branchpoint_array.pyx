
cdef class BranchpointArray:
    cdef readonly dict S_plus, S_minus, E_plus, E_minus, J_plus, J_minus
    cdef readonly int leftmost, rightmost, extend, end_extend, length, min_overhang
    cdef readonly tuple branchpoints
    cdef readonly float weight, threshold, cap_bonus, minimum_proportion
    cdef readonly np.ndarray depth, cov_plus, cov_minus
    cdef public OrderedArray bp_plus, bp_minus
    cdef readonly bint infer_starts, infer_ends, use_attributes
    def __init__(self, int leftmost, int rightmost, tuple reads, int extend, int end_extend, float minimum_proportion, float cap_bonus=5, int min_overhang=0, bint infer_starts=False, bint infer_ends=False, bint use_attributes=False):
        """Makes a collection of positions in the locus that act as dividing
        points for all nodes of the graph. Every donor (D) and acceptor (A) 
        site are retained as-is, but start (S) and end (E) are collapsed into
        a set of reference points."""
        cdef:
            BranchPoint BP, bp, lbp, rbp, lterm, rterm
            RNAseqMapping read
            float e_threshold, s_threshold, s_counter, e_counter, weight, wt, threshold_depth, s_weight, c_weight, e_weight
            str branchtype, junction_hash
            int pos, dist, bt, length, i, number_of_reads
            char strand
            list s_rank, e_rank, rank_sort, merged_plus, merged_minus, to_remove, donors, acceptors, gaps, gap_branchpoints
            (float, int, int, int) ranking
            (int, float) item
            (int, int) block
            dict lookup, end_counts, bp_dict 
            set donor_sites, acceptor_sites
            bint first_element, s_added
            tuple branchpoints
        
        bp_dict = {}
        self.J_plus = {}
        self.J_minus = {}
        self.extend = extend
        self.end_extend = end_extend
        self.minimum_proportion = minimum_proportion
        self.min_overhang = min_overhang
        self.cap_bonus = cap_bonus
        self.leftmost = leftmost
        self.rightmost = rightmost
        self.infer_starts = infer_starts
        self.infer_ends = infer_ends
        self.use_attributes = use_attributes
        self.threshold = 1 - self.minimum_proportion
        length = self.rightmost-self.leftmost
        self.depth_array(list(reads))
        
        # Populate the two BP arrays with D/A sites
        #TODO: Remove D/A sites below minimum_proportion of their position's coverage
        donors = self.bp_from_counter(Dp, 'D', 1)
        acceptors = self.bp_from_counter(Ap, 'A', 1)
        # Remove junctions from the lookup dicts if their donor and/or acceptor was filtered out
        donor_sites = set([bp.pos for bp in donors])
        acceptor_sites = set([bp.pos for bp in acceptors])
        for junction_hash in list(self.J_plus.keys()):
            block = literal_eval(junction_hash)
            if block[0] not in donor_sites or block[1] not in acceptor_sites:
                del self.J_plus[junction_hash]
        
        self.bp_plus = OrderedArray(donors + acceptors)
        
        donors = self.bp_from_counter(Dm, 'D', -1)
        acceptors = self.bp_from_counter(Am, 'A', -1)
        donor_sites = set([bp.pos for bp in donors])
        acceptor_sites = set([bp.pos for bp in acceptors])
        for junction_hash in list(self.J_minus.keys()):
            block = literal_eval(junction_hash)
            if block[1] not in donor_sites or block[0] not in acceptor_sites:
                del self.J_minus[junction_hash]
        
        self.bp_minus = OrderedArray(donors + acceptors)
        # Generate a priority queue for adding S/E branchpoints
        for strand in [1, -1]:
            if strand == -1:
                S, E, C = Sm, Em, Cm
            else:
                S, E, C = Sp, Ep, Cp

            s_threshold = sum(S.values()) * self.threshold
            e_threshold = sum(E.values()) * self.threshold
            
            s_rank = [
                (-S[pos], self.distance_from_edge(pos, strand, 'S'), 2, pos)
                for pos in sorted(S, reverse=True, key=S.get)
            ]
            e_rank = [
                (-E[pos], self.distance_from_edge(pos, strand, 'E'), 1, pos)
                for pos in sorted(E, reverse=True, key=E.get)
            ]
            s_counter = float(0)
            e_counter = float(0)
            rank_sort = sorted(s_rank + e_rank)
            for ranking in rank_sort:
                wt, dist, bt, pos = ranking
                weight = -wt
                if bt == 2:
                    if s_counter < s_threshold:
                        s_counter += weight
                        BP = StartPoint(strand, pos, weight, self.end_extend, C.get(pos, 0.0))
                        self.add_endpoint(BP, True)
                else:
                    if e_counter < e_threshold:
                        e_counter += weight
                        BP = EndPoint(strand, pos, weight, self.end_extend)
                        self.add_endpoint(BP, True)

        # After termination, remove all S/E branchpoints < minimum_proportion of their type
        end_counts = {branchtype:0.0 for branchtype in ['S','E','A','D','N']}
        for bp in self.bp_plus:
            if bp.branchtype == 'S':
                end_counts['S'] += bp.weight
            elif bp.branchtype == 'E':
                end_counts['E'] += bp.weight
            elif bp.branchtype in ['D','A']:
                if bp.weight > end_counts[bp.branchtype]:
                    end_counts[bp.branchtype] = bp.weight

        merged_plus = []
        first_element = True
        for bp in self.bp_plus:
            if bp.weight >= self.minimum_proportion * end_counts[bp.branchtype]:
                if bp.branchtype != 'S': # No extra restrictions for non-S branchpoints
                    merged_plus.append(bp)
                elif first_element:  # Is the leftmost branchpoint
                    merged_plus.append(bp)
                else:
                    merged_plus.append(bp)
                
                first_element = False       

        end_counts = {branchtype:0.0 for branchtype in ['S','E','A','D','N']}
        for bp in self.bp_minus:
            if bp.branchtype == 'S':
                end_counts['S'] += bp.weight
            elif bp.branchtype == 'E':
                end_counts['E'] += bp.weight
            elif bp.branchtype in ['D','A']:
                if bp.weight > end_counts[bp.branchtype]:
                    end_counts[bp.branchtype] = bp.weight

        merged_minus = []
        for  bp in self.bp_minus:
            if bp.weight >= self.minimum_proportion * end_counts[bp.branchtype]:
                if bp.branchtype != 'S': # No extra restrictions for non-S branchpoints
                    merged_minus.append(bp)
                else:
                    merged_minus.append(bp)
        
        # Allow the rightmost S- to be added regardless
        i = self.bp_minus.n - 1
        while i > 0:
            bp = self.bp_minus.ordered_array[i][0]
            if bp.weight < self.minimum_proportion * end_counts[bp.branchtype]: # BP will not be in bp_minus
                i -= 1
            elif bp.branchtype == 'S':
                if bp not in merged_minus:
                    merged_minus.append(bp)
                
                i = 0
            else: # A passing bp was found
                i = 0
        
        # Collapse all branchpoints into a single fixed tuple
        # Make a constant-time lookup to point to the appropriate branchpoint for any S/E read
        self.S_plus = {}
        self.S_minus = {}
        self.E_plus = {}
        self.E_minus = {}
        # Remove likely artifactual S/E branchpoints that are slightly overhanging from a heavier D/A element
        to_remove = []
        for i in range(len(merged_plus)):
            bp = merged_plus[i]
            if i > 0: # Check for an E < min_overhang right of a Donor site
                if bp.branchtype == 'E':
                    lbp = merged_plus[i-1]
                    if lbp.branchtype == 'D' and lbp.weight > bp.weight:
                        if bp.pos - lbp.pos < self.min_overhang:
                            to_remove.append(bp)
            
            if i + 1 < len(merged_plus):
                if bp.branchtype == 'S':
                    rbp = merged_plus[i+1]
                    if rbp.branchtype == 'A' and rbp.weight > bp.weight:
                        if rbp.pos - bp.pos < self.min_overhang:
                            to_remove.append(bp)
        
        merged_plus = [bp for bp in merged_plus if bp not in to_remove]

        to_remove = []
        for i in range(len(merged_minus)):
            bp = merged_minus[i]
            if i > 0: # Check for an E < min_overhang right of a Donor site
                if bp.branchtype == 'S':
                    lbp = merged_minus[i-1]
                    if lbp.branchtype == 'A' and lbp.weight > bp.weight:
                        if bp.pos - lbp.pos < self.min_overhang:
                            to_remove.append(bp)
            
            if i + 1 < len(merged_minus):
                if bp.branchtype == 'E':
                    rbp = merged_minus[i+1]
                    if rbp.branchtype == 'D' and rbp.weight > bp.weight:
                        if rbp.pos - bp.pos < self.min_overhang:
                            to_remove.append(bp)
        
        merged_minus = [bp for bp in merged_minus if bp not in to_remove]

        prohibited_positions = set()
        for bp in merged_plus:
            prohibited_positions.update(range(bp.pos-self.min_overhang, bp.pos+self.min_overhang+1))
        
        for bp in merged_minus:
            prohibited_positions.update(range(bp.pos-self.min_overhang, bp.pos+self.min_overhang+1))
        
        threshold_depth = np.max(self.depth)*self.minimum_proportion
        gaps = get_gaps(self.depth, self.extend, max(threshold_depth,1))
        gap_branchpoints = []
        for block in gaps:
            l = block[0]+self.leftmost
            r = block[1]+self.leftmost
            lbp = BranchPoint('N', 0, l, 0)
            rbp = BranchPoint('N', 0, r, 0)
            # Check if lbp and/or rbp should be added (>min_overhang from any existing branchpoint)
            if l not in prohibited_positions:
                gap_branchpoints.append(BranchPoint('N', 0, l, 0))
            
            if r not in prohibited_positions:
                gap_branchpoints.append(BranchPoint('N', 0, r, 0))


        branchpoints = tuple(sorted(merged_plus + merged_minus))
        # Evaluate whether terminal branchpoints must be added.
        # If no boundary Start/Endpoint exists, add a nonspecified one at leftmost/rightmost
        passes_threshold = np.where(self.depth >= threshold_depth)[0]
        if passes_threshold.shape[0] > 0:
            lbp = BranchPoint('N', 0, self.leftmost+passes_threshold[0], 0)
            rbp = BranchPoint('N', 0, self.leftmost+passes_threshold[-1], 0)
        else: # No positions pass the threshold
            lbp = BranchPoint('N', 0, self.leftmost, 0)
            rbp = BranchPoint('N', 0, self.rightmost, 0)
        
        if len(branchpoints) == 0:
            branchpoints = tuple([lbp, rbp])
        else:
            # Extend left border if S>/E< doesn't contain
            lterm = branchpoints[0]
            rterm = branchpoints[-1]
            add_lbp = True
            add_rbp = True
            if (lterm.branchtype == 'S' and lterm.strand == 1) or (lterm.branchtype == 'E' and lterm.strand == -1):
                if lbp.pos in lterm.span:
                    add_lbp = False

            # Extend right border if E>/S< doesn't contain
            if (rterm.branchtype == 'S' and rterm.strand == -1) or (rterm.branchtype == 'E' and rterm.strand == 1):
                if rbp.pos in rterm.span:
                    add_rbp = False
            
            if add_lbp or add_rbp:
                if add_lbp and add_rbp:
                    branchpoints = tuple([lbp] + list(branchpoints) + [rbp])
                elif add_lbp:
                    branchpoints = tuple([lbp] + list(branchpoints))
                else:
                    branchpoints = tuple(list(branchpoints) + [rbp])

        # Update 'N' branchpoints with their best inference
        self.branchpoints = self.infer_unknown_branchpoints(branchpoints)

        # Add an index for each branchpoint so a unique lookup can be used
        for i,bp in enumerate(self.branchpoints):
            bp.index = i
        
        # Add a lookup value for each S/E position to its corresponding branchpoint
        no_lookups = set(['D','A','N'])
        for bp in self.branchpoints:
            if bp.branchtype in no_lookups:
                continue
            
            if bp.branchtype == 'S':
                if bp.strand == 1:
                    bp_dict = Sp
                    lookup = self.S_plus
                else:
                    bp_dict = Sm
                    lookup = self.S_minus
            else:
                if bp.strand == 1:
                    bp_dict = Ep
                    lookup = self.E_plus
                else:
                    bp_dict = Em
                    lookup = self.E_minus
            
            i = bp.index
            span = bp.span
            # if i > 0:
            #     lb = self.branchpoints[bp.index - 1]
            #     span = range(max(lb.right, span.start), span.stop)
            
            # if i < len(self.branchpoints)-1:
            #     rb = self.branchpoints[bp.index + 1]
            #     span = range(span.start, min(rb.left+1, span.stop))
            
            for pos in bp_dict.keys():
                if pos in span:
                    if pos not in lookup:
                        lookup[pos] = i
                    else: # Check other branchpoints
                        otherBP = self.branchpoints[lookup[pos]]
                        # if abs(bp.pos-pos) < abs(otherBP.pos-pos): # Pick the closer one
                        if bp.weight > otherBP.weight: # Pick the heavier one
                            lookup[pos] = i
    
    cpdef void depth_array(self, list reads):
        """Returns a numpy array of coverage depth for a list of reads."""
        cdef:
            RNAseqMapping read
            int l, r, Sp, Cp, Ep, Dp, Ap, Sm, Cm, Em, Dm, Am, covp, covm, covn, covrow, pos
            float weight, s_weight, c_weight, e_weight
            (int, int) span, block
            np.ndarray coverage
            Py_ssize_t array_length
            str junction_hash
        
        # Populate a 13-column array:
        # S+  C+  E+  D+  A+  S-  C-  E-  D-  A-  cov+  cov-  cov.
        Sp, Cp, Ep, Dp, Ap, Sm, Cm, Em, Dm, Am, covp, covm, covn = range(13)
        array_length = self.rightmost - self.leftmost + 1 
        self.depth = np.zeros(shape=(13, array_length), dtype=np.float32)
        for read in reads:
            if self.use_attributes: # Check a read's attributes for different values of each type
                weight = read.weight
                s_weight = float(read.attributes.get('S.ppm', weight))
                e_weight = float(read.attributes.get('E.ppm', weight))
                c_weight = float(read.attributes.get('C.ppm', weight))
            else:
                weight = s_weight = e_weight = c_weight = read.weight
            
            if read.strand == 1:
                covrow = covp
                for span in read.junctions():
                    l = span[0] - self.leftmost
                    r = span[1] - self.leftmost
                    self.depth[Dp, l] += weight
                    self.depth[Ap, r] += weight
                    block = (l,r)
                    junction_hash = str(block)
                    self.J_plus[junction_hash] = self.J_plus.get(junction_hash, 0) + weight
                
                if read.s_tag:
                    pos = read.span[0] - self.leftmost
                    self.depth[Sp, pos] += s_weight
                    if read.capped:
                        self.depth[Cp, pos] += c_weight
                
                if read.e_tag:
                    pos = read.span[1] - self.leftmost
                    self.depth[Ep, pos] += e_weight
            elif read.strand == -1:
                covrow = covm
                for span in read.junctions():
                    l = span[0] - self.leftmost
                    r = span[1] - self.leftmost
                    self.depth[Am, l] += weight
                    self.depth[Dm, r] += weight
                    block = (l,r)
                    junction_hash = str(block)
                    self.J_minus[junction_hash] = self.J_minus.get(junction_hash, 0) + weight
                
                if read.e_tag:
                    pos = read.span[0] - self.leftmost
                    self.depth[Em, pos] += e_weight
                
                if read.s_tag:
                    pos = read.span[1] - self.leftmost
                    self.depth[Sm, pos] += s_weight
                    if read.capped:
                        self.depth[Cm, pos] += c_weight
            else: # The read has no features other than non-stranded coverage
                covrow = covn
            
            for span in read.ranges:
                l = span[0] - self.leftmost
                r = span[1] - self.leftmost
                self.depth[covrow, l:r] += weight
        
        # Multiply cap signal by 'cap_bonus'
        self.depth[Cp,:] *= self.cap_bonus
        self.depth[Cm,:] *= self.cap_bonus
        self.cov_plus = np.sum(self.depth[(covp,covn),:],axis=0)
        self.cov_minus = np.sum(self.depth[(covm,covn),:],axis=0)
    
    cpdef list bp_from_counter(self, dict counter, str branchtype, char strand):
        """Generates a list of Branchpoint objects
        of the given type and strand at the positions and weights
        given by a counter object. If the weight of the Branchpoint
        at pos fails to exceed minimum_proportion of coverage at
        that position it is filtered out."""
        cdef:
            int pos
            float weight
            BranchPoint bp
            list bp_list
        
        if branchtype == 'D': BP = DonorPoint
        elif branchtype == 'A': BP = AcceptorPoint
        elif branchtype == 'S': BP = StartPoint
        elif branchtype == 'E': BP = EndPoint
        else: BP = BranchPoint
        bp_list = []
        for pos in sorted(counter.keys()):
            weight = counter[pos]
            if weight >= self.minimum_proportion * self.depth[pos-self.leftmost]:
                bp = BP(strand, pos, weight)
                bp_list.append(bp)

        return bp_list

    cpdef int distance_from_edge(self, int pos, char strand, str branchtype):
        """ Given a position (int), strand (1,-1), and branchtype (S,E,D,A),
        returns the the distance to the most extreme edge of the locus."""
        if branchtype == 'S' or branchtype == 'D': # Consider upstream
            if strand == 1:
                return pos - self.leftmost
            else:
                return self.rightmost - pos
        else: # Consider downstream
            if strand == 1:
                return self.rightmost - pos
            else:
                return pos - self.leftmost
    
    cpdef void add_endpoint(self, BranchPoint BP, force=True):
        """Add an S/E BranchPoint to the bp dict of OrderedArrays.
        If at least one same-type branchpoint is in span, merge BP
        with the higher-priority branchpoint. If BP is flanked by two
        neighboring branchpoints, merge both BP and secondary into the first."""
        cdef OrderedArray array
        cdef BranchPoint left, right
        cdef int lp, rp, insert_site
        cdef bint l_in_range, r_in_range
        array = self.bp_plus if BP.strand == 1 else self.bp_minus        
        insert_site = array.probe(BP)
        lp = rp = -1
        l_in_range = r_in_range = False
        if insert_site > 0:
            left, lp = array.ordered_array[insert_site-1] # Get the left neighbor and its priority order
            if left.branchtype == BP.branchtype:
                if BP.pos in left.span: # BP could be merged with left
                    l_in_range = True
        
        if insert_site < array.n:
            right, rp = array.ordered_array[insert_site] # Get the right neighbor and its priority order
            if right.branchtype == BP.branchtype:
                if BP.pos in right.span: # BP could be merged with right
                    r_in_range = True
        
        if l_in_range:
            if r_in_range:
                if lp < rp: # left is prioritized over right
                    left.merge(BP)
                    left.merge(right)
                    array.delete(insert_site) # Remove right from the array
                else: # right is prioritized over left
                    right.merge(BP)
                    right.merge(left)
                    array.delete(insert_site-1) # Remove left from the array
            else: # only left can merge BP
                left.merge(BP)
        elif r_in_range: # only right can merge BP
            right.merge(BP)
        else: # no compatible neighbors
            if force:
                array.add(BP)
    
    cpdef tuple infer_unknown_branchpoints(self, tuple branchpoints):
        """Tries to use adjacent branchpoint features to infer the identity
        of gaps marked by 'N' (unknown) branchpoints. Returns a modified
        branchpoint tuple."""
        cdef tuple adjusted_branchpoints = branchpoints
        cdef BranchPoint FDP, FDM, LAP, LAM
        # If no inference is allowed, return with no changes
        if not self.infer_starts and not self.infer_ends:
            return adjusted_branchpoints
        
        # Iterate over branchpoints to identify where an inference is needed
        
        return adjusted_branchpoints

    cpdef void remove_branchpoint(self, str strand, int index):
        """Remove the BranchPoint at index in the strand bp dict.
        If the neighboring branchpoints are the same type and within range of each other,
        merge the smaller into the larger."""
        cdef OrderedArray array
        if strand == '+':
            array = self.bp_plus
        else:
            array = self.bp_minus
        
        array.delete(index)
        
        left, lp = array.ordered_array[index-1] # Get the left neighbor and its priority order
        if index < array.n:
            right, rp = array.ordered_array[index] # Get the right neighbor and its priority order
            if left.branchtype == right.branchtype: # Both neighbors were the same type
                if left.branchtype in ['S','E']: # Branchtype is terminal
                    if left.span[-1] >= right.span[0]: # Neighbor spans overlap
                        if lp < rp: # left is prioritized over right
                            left.merge(right)
                            array.delete(index) # Remove right from the array
                        else: # right is prioritized over left
                            right.merge(left)
                            array.delete(index-1) # Remove left from the array

