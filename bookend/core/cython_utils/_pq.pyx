#cython: language_level=3
cimport cython
import copy

'''
    Derived from Algorithms, 4th edition by Robert Sedgewick, Kevin Wayne
    See Coursera:  Algorithms, Part I - Week 4 - APIs and Elementary Implementations
    
    In binary heap a, a[0] is the root, which is the largest (maxPQ) or smallest (minPQ) key in the tree.
    The parent of a node at k is at position (k-1)/2
    
    Remember bitwise operators!
      x << y bitshift left by y (x * 2**y)
      x >> y bitshift right by y (x // 2**y) # // is floor division, do not use / in Python 3+
      x | y bitwise or ([1 if a or b else 0 for a,b in zip(bits(x),bits(y))])
      x & y bitwise and ([1 if a and b else 0 for a,b in zip(bits(x),bits(y))])
      ~ x bitflip (-x - 1)
      x ^ y bitwise exclusive or ([a if not b else ~a for a,b in zip(bits(x),bits(y))])
    It's possible to use the library heapq and add (priority, task) tuples to the heap.
    heapq documentation recommends (priority, entrycount, task) triples to break priority ties by insert order.
    It can be useful to store a discrete number of indexed items in a queue
    This way the queue never surpasses size V of the Graph
    See Algorithms 4e Chapter 2, p.320-322 for a full description of an indexed priority queue
'''

cdef class IndexMinPQ():
    """Stable implementation of an Indexed Minimum Priority Queue, 
    where each index is stored with its item and its insert order,
    and a lookup table allows traceback of each key's position on the queue.
    Can initialize on a list of (index, item) pairs."""
    cdef public list pq, qp
    cdef public int n, entry
    def __init__(self, N):
        self.pq = [] # Empty priority queue, stores triples e.g. ('<5--2 0.15>',0,5)
        self.qp = [None]*N # Reverse lookup dict that stores the current position of each index on the queue, e.g. 5:0
        self.n = 0 # Length of queue
        self.entry = 0 # Running total of items added
    
    def __len__(self): # Calling len(MinHeap) returns the queue length variable
        """Returns the length of the queue"""
        return self.n
    
    # _siftdown and _siftup are taken from the base Python library 'heapq',
    # with the exception that they now keep the reverse lookup dict updated when 
    # positions are swapped in the queue.
    cdef _siftdown(self, int startpos, int pos):
        cdef int parentpos
        cdef tuple newitem, parent
        newitem = self.pq[pos]
        # Follow the path to the root, moving parents down until
        # finding a place newitem fits.
        while pos > startpos:
            parentpos = (pos - 1) >> 1
            parent = self.pq[parentpos]
            if newitem < parent:
                self.pq[pos], self.qp[parent[-1]] = parent, pos # Update the queue and the reverse lookup dict
                pos = parentpos
                continue
            break
        
        self.pq[pos], self.qp[newitem[-1]] = newitem, pos # Update the queue and the reverse lookup dict
        
    cdef _siftup(self, int pos):
        cdef int endpos, startpos, childpos, rightpos
        cdef tuple newitem
        endpos = len(self.pq)
        startpos = pos
        newitem = self.pq[pos]
        # Bubble up the smaller child until hitting a leaf.
        childpos = 2*pos + 1    # leftmost child position
        while childpos < endpos:
            # Set childpos to index of smaller child.
            rightpos = childpos + 1
            if rightpos < endpos and not self.pq[childpos] < self.pq[rightpos]:
                childpos = rightpos
            # Move the smaller child up.
            self.pq[pos], self.qp[self.pq[childpos][-1]] = self.pq[childpos], pos # Update the queue and the reverse lookup dict
            
            pos = childpos
            childpos = 2*pos + 1
        # The leaf at pos is empty now.  Put newitem there, and bubble it up
        # to its final resting place (by sifting its parents down).
        self.pq[pos], self.qp[newitem[-1]] = newitem, pos # Update the queue and the reverse lookup dict
        self._siftdown(startpos, pos)
    
    cdef heappush(self, tuple q_item):
        """Push q_item onto heap, maintaining the heap invariant."""
        self.pq.append(q_item)
        self.qp[q_item[-1]] = len(self.pq)-1
        self._siftdown(0, len(self.pq)-1)

    cdef tuple heappop(self):
        """Pop the smallest item off the heap, maintaining the heap invariant."""
        cdef tuple lastelt, returnitem
        lastelt = self.pq.pop()    # raises appropriate IndexError if heap is empty
        if self.pq:
            returnitem = self.pq[0]
            self.pq[0] = lastelt
            self.qp[lastelt[-1]] = None
            self._siftup(0)
            return returnitem
        
        return lastelt
    
    def insert(self, index, item):
        """Add a single indexed item to an existing priority queue"""
        if self.qp[index] is not None: # This index is already stored on the queue
            current_loc = self.qp[index]
            current_value = self.pq[current_loc][0] # Get the existing value for this key
            if current_value <= item: # The priority of the old data is higher, do not add
                return
            
            # print("**UPDATE**")
            q_item = (item, self.entry, index) # Make a tuple that can be sifted through the heap
            self.pq[current_loc] = q_item # Do an in-place substitution of this index's q_item
            self._siftdown(0, current_loc) # Swim the updated item to its place in the queue
        else: # This key doesn't exist yet. Add the item to the queue
            # print("**NEW**")
            self.entry += 1
            self.n += 1
            q_item = (item, self.entry, index) # Make a tuple that can be sifted through the heap
            self.heappush(q_item) # Push item onto heap and get its position in the queue
        
        # assert self.isInvariant()
    
    def contains(self, index):
        """Check whether a certain index is on the queue"""
        if self.qp[index] is None:
            return False
        
        return True
    
    def pop(self, withindex=False):
        """Return and remove the smallest key"""
        cdef int index
        nextitem = self.pq[0] # The first element is the min
        if self.n == 1:
            self.pq = []
            self.n = 0
        else:
            nextitem = self.heappop() # Remove the last item to position 0 of the heap
            self.n += -1 # Reduce the size indicator of the heap
        
        item, index = nextitem[0], nextitem[-1]
        self.qp[index] = None # Remove this index from the lookup dict
        if withindex:
            return (index, item)
        
        return item
    
    def min(self):
        """Peek at the top item in the queue without removing"""
        return self.pq[0][0]
    
    def minIndex(self):
        """Peek at the top item's value without removing"""
        return self.pq[0][-1]
    
    def isEmpty(self):
        """Boolean: is the queue empty?"""
        if self.n == 0:
            return True
        else:
            return False
    
    def dump(self,withindex=False):
        """Return all keys in the heap in ascending order"""
        output = []
        while self.n > 0:
            output += [self.pop(withindex)]
        
        return output
    
    def isInvariant(self):
        """Perform a check that the heap is in fact invariant."""
        current_heap = copy.deepcopy(self)
        dequeued = current_heap.dump()
        # print([a for a,b,c in sorted(self.pq)])
        # print(dequeued)
        return [a for a,b,c in sorted(self.pq)] == dequeued
        
