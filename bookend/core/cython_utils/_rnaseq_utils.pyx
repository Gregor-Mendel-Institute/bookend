#cython: language_level=3
cimport cython
import array
from cpython cimport array
import numpy as np
cimport numpy as np
import json
from cython_utils import _fasta_utils as fu
from collections import namedtuple, Counter
from ast import literal_eval
import copy
ctypedef unsigned char uint8
ctypedef np.float32_t float32

##############################################################################################################
##############################################################################################################
# Object for defining a collection of RNA sequencing reads #
##############################################################################################################
##############################################################################################################
config_defaults = {
    'source':'',
    's_tag':False,
    'e_tag':False,
    'capped':False,
    'stranded':False,
    'start_seq':'ACGGG',
    'end_seq':'RRRRRRRRRRRRRRR',
    'minlen':20,
    'mismatch_rate':0.2,
    'sj_shift':2
}

cdef class RNAseqDataset:
    cdef public list read_list, chrom_array, source_array, chrom_lengths
    cdef public int chrom_index, source_index
    cdef public dict chrom_dict, source_dict
    cdef readonly dict config, genome
    cdef readonly bint s_tag, e_tag, capped, stranded, ignore_ends
    cdef readonly str start_seq, end_seq
    cdef readonly int minlen, sj_shift
    cdef readonly float mismatch_rate
    cdef readonly array.array start_array, end_array

    def __init__(self, chrom_array=None, source_array=None, chrom_lengths=None, genome_fasta=None, config=config_defaults):
        """Container for RNAseqMapping objects. Stores a reference dictionary for all
        chromosome names and sample names. Contains methods for parsing
        a variety of files into a collection of read objects."""
        self.read_list = []
        self.config = config
        self.s_tag = self.config['s_tag']
        self.e_tag = self.config['e_tag']
        self.capped = self.config['capped']
        self.stranded = self.config['stranded']
        self.start_seq = self.config['start_seq']
        self.end_seq = self.config['end_seq']
        self.minlen = self.config['minlen']
        self.mismatch_rate = self.config['mismatch_rate']
        self.sj_shift = self.config['sj_shift']
        self.start_array = fu.nuc_to_int(self.start_seq)
        self.end_array = fu.nuc_to_int(self.end_seq)
        if genome_fasta is not None:
            self.genome = fu.import_genome(genome_fasta)
            self.chrom_array = sorted(self.genome.keys())
            self.chrom_index = len(self.chrom_array)
            self.chrom_dict = dict(zip(self.chrom_array, range(self.chrom_index)))
            self.chrom_lengths = [len(self.genome[chrom]) for chrom in self.chrom_array]
        else:
            self.genome = {}
            self.chrom_lengths = chrom_lengths
            self.chrom_dict = {}
            self.chrom_index = 0
            self.chrom_array = []
            if chrom_array is not None:
                for c in chrom_array:
                    self.add_chrom(c)
        
        self.source_dict = {}
        self.source_index = 0
        self.source_array = []
        if source_array is not None:
            for s in source_array:
                self.add_source(s)
    
    cpdef add_source(self, source_string):
        if source_string not in self.source_dict:
            self.source_array.append(source_string)
            self.source_dict[source_string] = self.source_index
            self.source_index += 1
    
    cpdef add_chrom(self, chrom_string):
        if chrom_string not in self.chrom_dict:
            self.chrom_array.append(chrom_string)
            self.chrom_dict[chrom_string] = self.chrom_index
            self.chrom_index += 1
    
    cpdef add_read_from_BED(self, bed_line, source_string=None, s_tag=False, e_tag=False, capped=False, gaps_are_junctions=False):
        cdef RNAseqMapping new_read
        input_data = parse_BED_line(bed_line, self.chrom_dict, self.source_dict, source_string, s_tag, e_tag, capped, gaps_are_junctions)
        if type(input_data.chrom) is str: # Chromosome wasn't in chrom_dict
            chrom_string = input_data.chrom
            input_data = input_data._replace(chrom=self.chrom_index)
            self.add_chrom(chrom_string)
        
        if type(input_data.source) is str: # Source wasn't in source_dict
            source_string = input_data.source
            input_data = input_data._replace(source=self.source_index)
            self.add_source(source_string)
        
        new_read = RNAseqMapping(input_data)
        self.read_list.append(new_read)
    
    cpdef add_read(self, input_data):
        cdef RNAseqMapping new_read
        if type(input_data.chrom) is str: # Chromosome wasn't in chrom_dict
            chrom_string = input_data.chrom
            input_data = input_data._replace(chrom=self.chrom_index)
            self.add_chrom(chrom_string)
        
        if type(input_data.source) is str: # Source wasn't in source_dict
            source_string = input_data.source
            input_data = input_data._replace(source=self.source_index)
            self.add_source(source_string)
        
        new_read = RNAseqMapping(input_data)
        self.read_list.append(new_read)
    
    cpdef add_read_from_ELR(self, elr_line):
        cdef RNAseqMapping new_read
        new_read = elr_to_readobject(elr_line)
        self.read_list.append(new_read)
    
    cpdef add_read_from_BAM(self, bam_lines, bint ignore_ends=False, bint secondary=False):
        cdef list new_read_list
        if type(bam_lines) is not list:
            bam_lines = [bam_lines]
        
        new_read_list = generate_read_from_bam(self, bam_lines, ignore_ends, secondary)
        self.read_list += new_read_list
    
    cpdef add_read_from_GTF(self, gtf_lines, bint ignore_ends=False):
        cdef RNAseqMapping new_read
        new_read = generate_read_from_gtf(self, gtf_lines)
        self.read_list += new_read
    
    cpdef add_read_from_GFF(self, gff_lines, bint ignore_ends=False):
        cdef RNAseqMapping new_read
        new_read = generate_read_from_gff(self, gff_lines)
        self.read_list += new_read

    cpdef pop_read(self, read_format='elr', as_string=True):
        """Remove the last read added to the stack and write it in 'format'.
        """
        if read_format.lower() == 'elr':
            return(self.read_list.pop().write_as_elr(as_string))
        elif read_format.lower() == 'bed':
            return(self.read_list.pop().write_as_bed(self.chrom_array, self.source_array, as_string))
        elif read_format.lower() == 'gtf':
            return(self.read_list.pop().write_as_gtf(self.chrom_array, 'bed'))
    
    cpdef dump_header(self):
        """Returns an array of strings that describe chrom_dict and source_dict of the Dataset."""
        header_list = []
        for i,c in enumerate(self.chrom_array):
            header_list += ['#C {} {}'.format(i, c)]
        
        for i,s in enumerate(self.source_array):
            header_list += ['#S {} {}'.format(i, s)]
        
        return header_list

##############################################################################################################
##############################################################################################################
# Object for defining an RNA sequencing read #
##############################################################################################################
##############################################################################################################
ELdata = namedtuple('ELdata', 'chrom source strand ranges splice s_tag e_tag capped weight')
cdef class RNAseqMapping:
    cdef public int chrom, source, strand
    cdef public list ranges, splice
    cdef public bint s_tag, e_tag, capped, complete
    cdef public dict attributes
    cdef public (int, int) span
    cdef public float weight, coverage
    def __init__(self, input_data, attributes = None):
        """Initializes a Read Object given a tuple of input data.
        Requires a chromosome, strand, source, weight, a sorted tuple of
        exon ranges and an array of booleans indicating which gaps between exons are splice junctions."""
        self.chrom, self.source = int(input_data.chrom), int(input_data.source)
        self.strand = input_data.strand
        self.ranges = input_data.ranges
        self.splice = input_data.splice
        if self.strand == 0: # Terminal tag information is meaningless for nonstranded reads
            self.s_tag = self.e_tag = self.capped = False
        else:
            self.s_tag, self.e_tag, self.capped = input_data.s_tag, input_data.e_tag, input_data.capped
        
        self.weight = float(input_data.weight)
        self.span = (self.left(), self.right())
        self.complete = False
        if self.s_tag and self.e_tag and False not in self.splice:
            self.complete = True
        
        if attributes is not None:
            self.attributes = attributes
        else:
            self.attributes = {}
    
    def __eq__(self, other): return self.span == other.span
    def __gt__(self, other): return self.span > other.span
    def __ge__(self, other): return self.span >= other.span
    def __lt__(self, other): return self.span < other.span
    def __le__(self, other): return self.span <= other.span
    def __ne__(self, other): return self.span != other.span
    
    def __repr__(self):
        """Represents the read object with an ASCII character string:
            >>, <<, || represent plus, minus, and unstranded
            Ranges are connected by ^ (splice) or . (gap)"""
        if self.strand == 1:
            strandchar = '>'
        elif self.strand == -1:
            strandchar = '<'
        else:
            strandchar = '|'
        
        gapchar = ['^' if i else '_' for i in self.splice] + [strandchar]
        rangechar = ['{}-{}'.format(a,b) for a,b in self.ranges]
        return(''.join([a+b for a,b in zip(rangechar,gapchar)]))
    
    def __len__(self):
        return self.right() - self.left()
    
    cpdef int left(self):
        return self.ranges[0][0]
    
    cpdef int right(self):
        return self.ranges[-1][-1]
    
    cpdef int get_length(self):
        """Returns the number of nucleotides covered by all blocks of the object."""
        cdef int length = 0
        cdef (int,int) exon
        for exon in self.ranges:
            length += exon[1] - exon[0]
        
        return length

    def gaps(self):
        """Returns an array of 0-indexed (start, end) tuples of gaps between ranges"""
        if len(self.ranges) == 1:
            return []
        
        return [(self.ranges[i][-1], self.ranges[i+1][0]) for i in range(len(self.ranges)-1)]
    
    def junctions(self):
        """Returns an array of 0-indexed (start, end) tuples of intron locations"""
        j_array = []
        for i,j in enumerate(self.splice):
            if j:
                j_array += [(self.ranges[i][-1], self.ranges[i+1][0])]
        
        return j_array
    
    cpdef bint overlaps(self, RNAseqMapping other):
        """Returns a boolean if the mapping range of self overlaps other."""
        cdef int l1, r1, l2, r2
        l1, r1 = self.span
        l2, r2 = other.span
        if r1 < l2 or l1 > r2:
            return False
        
        return True
    
    cpdef (int, int) overlap_range(self, RNAseqMapping other):
        """Returns 0-indexed open range of overlap between self and other"""
        edges = sorted([self.span, other.span])
        return (edges[1][0], edges[0][1])
    
    cpdef bint ends_clash(self, RNAseqMapping other):
        """Returns a boolean of whether the combination of end tags between
        self and other can be substrings of the same object."""
        cdef (int, int) strands = tuple(sorted([self.strand, other.strand]))
        if strands[0] == 1:
            assert strands[1] != -1
            if self.s_tag and self.span[0] > other.span[0]: return True # Other is left of self's S tag
            if other.s_tag and other.span[0] > self.span[0]: return True # Self is left of other's S tag
            if self.e_tag and self.span[1] < other.span[1]: return True # Other is right of self's E tag
            if other.e_tag and other.span[1] < self.span[1]: return True # Self is right of other's E tag
        elif strands[0] == -1:
            assert strands[1] != 1
            if self.s_tag and self.span[1] < other.span[1]: return True # Other is right of self's S tag
            if other.s_tag and other.span[1] < self.span[1]: return True # Self is right of other's S tag
            if self.e_tag and self.span[0] > other.span[0]: return True # Other is left of self's E tag
            if other.e_tag and other.span[0] > self.span[0]: return True # Self is left of other's E tag
        
        return False
    
    cpdef bint is_compatible(self, RNAseqMapping other):
        """Self and other contain no attributes that demonstrate they could
        not be subsequences of a common longer molecule."""
        cdef (int, int) overlap
        if self.chrom != other.chrom:
            return False
        
        if self.source != other.source:
            return False
        
        if self.strand != 0 and other.strand != 0 and self.strand != other.strand:
            return False
        
        # If self or other contain terminal tags, the tags must be compatible
        if self.ends_clash(other):
            return False
        
        if not self.overlaps(other):
            return True # No incompatibilities were found
        
        # If the two reads share a chrom and strand and overlap,
        # check the overlapping range for identical splice architecture
        overlap = self.overlap_range(other)
        j1 = [j for j in self.junctions() if j[0] > overlap[-1] and j[-1] < overlap[0]]
        j2 = [j for j in other.junctions() if j[0] > overlap[-1] and j[-1] < overlap[0]]
        if j1 == j2:
            return True
        
        return False
        
    cpdef bint is_identical(self, RNAseqMapping other):
        """Boolean of whether two read objects share all nontrivial attributes"""
        if self.source != other.source: return False
        if self.chrom != other.chrom: return False
        if self.strand != other.strand: return False
        if self.ranges != other.ranges: return False
        if self.splice != other.splice: return False
        if self.s_tag != other.s_tag: return False
        if self.e_tag != other.e_tag: return False
        if self.capped != other.capped: return False
        return True
    
    cpdef bint merge(self, RNAseqMapping other):
        """Combines with another Mapping object. Must be compatible."""
        if not self.is_compatible(other):
            return False
        
        # Unify the strand information of the two objects
        self.s_tag = self.s_tag or other.s_tag
        self.e_tag = self.e_tag or other.e_tag
        self.capped = self.capped or other.capped
        self.weight = self.weight + other.weight
        if self.strand == 0:
            self.strand = other.strand
        
        if self.overlaps(other): # The two ranges overlap to some degree
            junctions = self.junctions() + other.junctions()
            self.ranges = collapse_blocks(self.ranges + other.ranges)
            new_gaps = self.gaps()
            self.splice = [gap in junctions for gap in new_gaps]
        elif self < other: # self is strictly left of other
            self.ranges = self.ranges + other.ranges
            self.splice = self.splice + [False] + other.splice # Join the two ranges with a gap
        else: # self is strictly right of other
            self.ranges = other.ranges + self.ranges
            self.splice = other.splice + [False] + self.splice # Join the two ranges with a gap
        
        self.span = (self.left(), self.right())
        return True
    
    def get_node_labels(self):
        """Returns a string with one label for each edge of each range in self.ranges."""
        if self.strand == 1:
            gapchar = 'DA'
            if self.s_tag:
                if self.capped:
                    startchar = 'C'
                else:
                    startchar = 'S'
            else:
                startchar = '.'
            
            if self.e_tag:
                endchar = 'E'
            else:
                endchar = '.'
        elif self.strand == -1:
            gapchar = 'AD'
            if self.s_tag:
                if self.capped:
                    endchar = 'C'
                else:
                    endchar = 'S'
            else:
                endchar = '.'
            
            if self.e_tag:
                startchar = 'E'
            else:
                startchar = '.'
        else:
            gapchar = '..'
            startchar = endchar = '.'
        
        return(''.join([startchar]+[gapchar if i else '..' for i in self.splice]+[endchar]))
    
    def write_as_elr(self, as_string=True):
        """Returns a string that represents the ReadObject
        in the end-labeled read (ELR) format"""
        elr_strand = '.'
        if self.strand == 1:
            elr_strand = '+'
        elif self.strand == -1:
            elr_strand = '-'
        
        block_ends = flatten(self.ranges)
        lengths = [block_ends[i]-block_ends[i-1] for i in range(1,len(block_ends))]
        labels = self.get_node_labels()
        EL_CIGAR = ''.join([str(a)+str(b) for a,b in zip(labels,lengths+[''])])
        read_len = self.right() - self.left()
        elr_line = [self.chrom, self.left(), read_len, elr_strand, EL_CIGAR, self.source, round(self.weight,2)]
        if as_string:
            return '\t'.join([str(i) for i in elr_line])
        else:
            return elr_line
     
    def write_as_bed(self, chrom_array, source_array, as_string=True, score_column='weight'):
        """Returns a string that represents the ReadObject
        in a 15-column BED format"""
        labels = self.get_node_labels()
        bed_strand = '.'
        if self.strand == 1:
            bed_strand = '+'
        elif self.strand == -1:
            bed_strand = '-'
        
        l = labels[0]
        r = labels[-1]
        ends = l+r
        if ends in ['SE','ES']:
            rgb = bed_colors['SE']
        elif 'C' in ends:
            rgb = bed_colors['C']
        elif 'S' in ends:
            rgb = bed_colors['S']
        elif 'E' in ends:
            rgb = bed_colors['E']
        else:
            rgb = bed_colors['U']
        
        name = '.'
        if score_column == 'weight':
            score = round(self.weight,2)
        elif score_column == 'coverage':
            score = round(self.coverage,2)
        else:
            score = '.'
        
        chromStart, blockStarts, blockSizes = explode_block_ranges(self.ranges)
        bed_line = [
            chrom_array[self.chrom], chromStart, self.right(),
            name, score, bed_strand, chromStart, self.right(), rgb,
            len(self.ranges),
            ','.join([str(i-1) for i in blockSizes]),
            ','.join([str(i) for i in blockStarts]),
            round(self.weight,2), source_array[self.source], labels
        ]
        if as_string:
            return '\t'.join([str(i) for i in bed_line])
        else:
            return bed_line
    
    def write_as_gtf(self, list chrom_array, str source):
        """Returns a string representation of the ReadObject as a set
        of 'transcript' and 'exon' lines in a GTF, with attributes
        passed in the 'transcript' line."""
        #TODO: Implement export as GTF
        cdef str gtf_txt = ''
        cdef str chrom = chrom_array[self.chrom]
        cdef str strand = {1:'+',-1:'-',0:'.'}[self.strand]
        if self.attributes is None: self.attributes = {}
        frame = '.'
        cdef str score = str(round(self.coverage,2))
        gtf_txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, source, 'transcript', self.left()+1, self.right(), score, strand, frame, self.print_attributes('gtf',True))
        for left,right in self.ranges:
           gtf_txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, source, 'exon', left+1, right, score, strand, frame, self.print_attributes('gtf',False))
        
        return gtf_txt[:-1]

    def print_attributes(self, format='gtf', full=True):
        """Outputs a text string of the key:value pairs in self.attributes,
        formatted as GTF or GTF."""
        if full:
            attributes = self.attributes
        else:
            attributes = {k:self.attributes[k] for k in ('gene_id','transcript_id')}
        
        if format == 'gtf':
            return ' '.join(['{} "{}";'.format(k,v) for k,v in attributes.items()])
        elif format == 'gff':
            return ' '.join(['{}={};'.format(k,v) for k,v in attributes.items()])[:-1]

    def write_as_json(self, chrom_array=None, source_array=None):
        """Returns a string that represents the ReadObject
        in a JSON format"""
        if chrom_array is None:
            chrom = self.chrom
        else:
            chrom = chrom_array[self.chrom]
        
        if source_array is None:
            source = self.source
        else:
            source = source_array[self.source]
        
        json_line = f'{{ \
        "chrom": {json.dumps(chrom)}, \
        "strand": {json.dumps(self.strand)} \
        "span": {json.dumps(self.span)}\
        "ranges": {json.dumps(self.ranges)}\
        "splice": {json.dumps(self.splice)}\
        "s_tag": {json.dumps(self.s_tag)}\
        "e_tag": {json.dumps(self.e_tag)}\
        "capped": {json.dumps(self.capped)}\
        "weight": {json.dumps(self.weight)}\
        "source": {json.dumps(source)}\
        }};'
        
        return json_line

bed_colors = {
    'U':'128,130,133',
    'C':'0,212,145',
    'S':'28,117,188',
    'E':'190,30,45',
    'SE':'109,74,116'
}

def array_to_blocks(list arr):
    """Breaks an array into a collections of (start,end) INCLUSIVE doubles
    that define the coordinates of all contiguous blocks in the array"""
    cdef list clean_array
    cdef int block_start, block_end, i
    clean_array = sorted(list(set([i for i in arr if type(i) is int])))
    if len(clean_array) == 0:
        raise StopIteration
    
    block_start = clean_array[0]
    block_end = clean_array[0]
    if len(arr) == 1:
        yield (block_start, block_end)
        raise StopIteration
    
    for i in clean_array[1:]:
        if i-1 == block_end:
            block_end = i
        else:
            yield (block_start, block_end)
            block_start = i
            block_end = i
    
    yield (block_start, block_end)

cpdef tuple overlap_type(tuple range_a, tuple range_b):
    """Returns a list of relationships
    between the edges of range_a and range_b"""
    cdef int LL, LR, RL, RR
    LL = range_a[0] >= range_b[0]
    LR = range_a[0] >= range_b[1]
    RL = range_a[1] >= range_b[0]
    RR = range_a[1] >= range_b[1]
    cdef tuple overlap_pattern = (LL, LR, RL, RR)
    return overlap_pattern

cpdef list flatten(list list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

cpdef list collapse_blocks(list list_of_doubles):
    """Iterates over a sorted list of 0-indexed open (start,end) doubles
    and merges those that overlap to return a set of nonoverlapping blocks"""
    cdef list new_blocks
    cdef (int, int) current_block, previous_block
    cdef bint previous_block_exists
    cdef int overlap, merged_left, merged_right
    list_of_doubles.sort()
    new_blocks = []
    previous_block_exists = False
    for current_block in list_of_doubles:
        if previous_block_exists:
            # Get the overlap type of the two adjacent blocks
            overlap = sum(overlap_type(current_block, previous_block))
            if overlap > 0 and overlap < 4:
                # The two blocks do overlap to some degree
                merged_left = min(previous_block[0],current_block[0])
                merged_right = max(previous_block[-1],current_block[-1])
                current_block = (merged_left, merged_right)
            else:
                new_blocks += [previous_block]
        
        # Update the last block
        previous_block = current_block
        if not previous_block_exists:
            previous_block_exists = True
    
    # Resolve the last block
    new_blocks += [previous_block]
    return new_blocks

cdef list get_block_ranges(str chromStart, str blockStarts, str blockSizes):
    """Reorganizes the BED12 block format into a list of 0-indexed open (start,end) doubles."""
    cdef list sizes, starts, block_ranges
    cdef int c, i, number_of_blocks, l, r, start, size
    cdef str s
    sizes = [int(s) for s in blockSizes.split(',')]
    starts = [int(s) for s in blockStarts.split(',')]
    c = int(chromStart)
    number_of_blocks = len(sizes)
    block_ranges =  [()] * number_of_blocks
    for i in range(number_of_blocks):
        start = starts[i]
        size = sizes[i]
        l = c+start
        r = l+size
        block_ranges[i] = (l, r)
    
    return block_ranges

cdef explode_block_ranges(block_ranges):
    """Converts a list of block ranges generated by get_block_ranges()
    back into a set of three variables: chromStart, blockStarts, blockSizes"""
    chromStart = block_ranges[0][0] # Leftmost position
    blockStarts = [a-chromStart for a,b in block_ranges]
    blockSizes = [b-a+1 for a,b in block_ranges]
    return chromStart, blockStarts, blockSizes

cdef parse_ELR_line(str elr_line):
    """Parses one line of an ELR file into an ELdata namedtuple.
    Example:
      0 6787 1950 - E282A87D294A113D86A586D90A91D48A129D144S 0 1.0
    """
    cdef str chrom_num, chromStart, read_len, elr_strand, EL_CIGAR, source_num, weight_string
    cdef int chrom, source, startpos, strand, position, i, f, a, b, rightside, l_index
    cdef float weight
    cdef list label_indices, feature_lengths, ranges, splice
    cdef bint gap
    chrom_num, chromStart, read_len, elr_strand, EL_CIGAR, source_num, weight_string = elr_line.rstrip().split('\t')
    if elr_strand == '+':
        strand = 1
    elif elr_strand == '-':
        strand = -1
    else:
        strand = 0
    
    chrom = int(chrom_num)
    source = int(source_num)
    weight = round(float(weight_string),2)
    startpos = int(chromStart)
    label_indices = [i for i,character in enumerate(EL_CIGAR) if not character.isdigit()]
    feature_lengths = [int(EL_CIGAR[a+1:b]) for a,b in zip(label_indices[:-1],label_indices[1:])]
    ranges = []
    splice = []
    gap = False
    position = startpos
    for i,f in enumerate(feature_lengths):
        if not gap:
            rightside = position + f
            ranges += [(position, rightside)]
            position = rightside
        else:
            position += f
            l_index = label_indices[i]
            if EL_CIGAR[l_index] in ['A','D']:
                splice += [True]
            else:
                splice += [False]
        
        gap = not gap
    
    s_tag = False
    e_tag = False
    capped = False
    if EL_CIGAR[0] == 'C':
        s_tag = True
        capped = True
    elif EL_CIGAR[0] == 'S':
        s_tag = True
    elif EL_CIGAR[-1] == 'C':
        s_tag = True
        capped = True
    elif EL_CIGAR[-1] == 'S':
        s_tag = True
    
    if EL_CIGAR[0] == 'E':
        e_tag = True
    elif EL_CIGAR[-1] == 'E':
        e_tag = True
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)

cdef parse_BED_line(bed_line, chrom_dict, source_dict, source_string=None, s_tag=False, e_tag=False, capped=False, gaps_are_junctions=False):
    """Parses one line of a 12- or 15-column BED file into an ELdata namedtuple.
    Examples:
      Ath_chr1 6787 8737 AT1G01020.6 . - 6787 8737 0,0,0 6 282,294,86,90,48,144 0,369,776,1448,1629,1806
      Ath_chr1 6787 8737 . . - 0 0 109,74,116 6 282,294,86,90,48,144 0,369,776,1448,1629,1806 1.0 TAIR10.40 EADADADADADS
    """
    bed_elements = bed_line.rstrip().split('\t')
    if len(bed_elements) == 15:
        chrom_string, chromStart, end, readname, score, bed_strand, mmnum, mmorder, rgb, blocknum, blockSizes, blockStarts, weight, source_string, label = bed_elements
    elif len(bed_elements) == 12:
        chrom_string, chromStart, end, readname, score, bed_strand, mmnum, mmorder, rgb, blocknum, blockSizes, blockStarts = bed_elements
        try:
            weight = float(score)
        except:
            weight = float(1)
        
        label = None
    else:
        return
    
    if bed_strand == '+':
        strand = 1
    elif bed_strand == '-':
        strand = -1
    else:
        strand = 0
    
    ranges = get_block_ranges(chromStart, blockStarts, blockSizes)
    if label is not None: # Determine what kind of end labels exist based on the label
        s_tag = e_tag = capped = False
        if strand == 1:
            if label[0] == 'C':
                s_tag = capped = True
            elif label[0] == 'S':
                s_tag = True
            
            if label[-1] == 'E':
                e_tag = True
            
            splice = [True if i=='D' else False for i in label[1:-1:2]]
        elif strand == -1:
            if label[-1] == 'C':
                s_tag = capped = True
            elif label[-1] == 'S':
                s_tag = True
            
            if label[0] == 'E':
                e_tag = True
            
            splice = [True if i=='A' else False for i in label[1:-1:2]]
        else:
            splice = [False]*(len(ranges)-1)
    else:
        if gaps_are_junctions:
            splice = [True]*(len(ranges)-1)
        else:
            splice = [False]*(len(ranges)-1)
    
    if chrom_string in chrom_dict:
        chrom = chrom_dict[chrom_string]
    else:
        chrom = chrom_string
    
    if source_string is not None:
        if source_string in source_dict:
            source = source_dict[source_string]
        else:
            source = source_string
    else:
        source = None
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)

cpdef RNAseqMapping elr_to_readobject(elr_line):
    """Converts an ELR line to an RNAseqMapping object
    """
    cdef input_data = parse_ELR_line(elr_line)
    output_object = RNAseqMapping(input_data)
    return output_object

##############################################################################################################
# GTF/GFF file parsing #
##############################################################################################################
gtf_defaults = {
    'parent_types':set(['transcript']),
    'parent_key_transcript':'transcript_id',
    'parent_key_gene':'gene_id',
    'child_types':set(['exon']),
    'child_key_transcript':'transcript_id',
    'child_key_gene':'gene_id'
}
gff_defaults = {
    'parent_types':set(['mRNA','transcript']),
    'parent_key_transcript':'ID',
    'parent_key_gene':'Parent',
    'child_types':set(['exon']),
    'child_key_transcript':'Parent',
    'child_key_gene':'gene_id'
}

cpdef dict parse_attributes(str attr_string, str attr_format):
    """Converts a GTF or GFF3 formatted string to """
    cdef dict attr_dict = {}
    if attr_format == 'GTF':
        attr_dict = {k:v for k,v in [attr.rstrip('";').split(' "') for attr in attr_string.split('; ')]}
    elif attr_format == 'GFF':
        attr_dict = {k:v for k,v in [attr.rstrip(';').split('=') for attr in attr_string.split(';')]}
    
    return attr_dict

cdef class AnnotationObject:
    cdef public dict attributes
    cdef public bint keep, parent
    cdef public str format, gene_id, transcript_id, chrom, source, anno_type
    cdef public int strand
    cdef public (int, int) span
    #cdef public (str, str, str, str, str, str, str, str, str) fields
    cdef public tuple fields
    def __init__(self, anno_string, format, config_dict):
        """Generates an intermediate object from a single line of a GTF/GFF file.
        This will not completely represent """
        cdef:
            set child_types, parent_types
            str gene_id_key, transcript_id_key
        
        self.format = format
        self.keep = False
        self.fields = tuple(anno_string.split('\t')[0:9])
        self.anno_type = self.fields[2]
        self.parent = False
        child_types = config_dict['child_types']
        parent_types = config_dict['parent_types']
        if self.anno_type in child_types:
            self.keep = True
            self.parent = False
        elif self.anno_type in parent_types:
            self.keep = True
            self.parent = True
        
        if self.keep:
            self.chrom = self.fields[0]
            self.source = self.fields[1]
            self.span = (int(self.fields[3])-1, int(self.fields[4]))
            self.attributes = parse_attributes(self.fields[8], self.format)
            if self.parent:
                gene_id_key=config_dict['parent_key_gene']
                transcript_id_key=config_dict['parent_key_transcript']
            else:
                gene_id_key=config_dict['child_key_gene']
                transcript_id_key=config_dict['child_key_transcript']
            
            self.gene_id = self.attributes.get(gene_id_key,'')
            self.transcript_id = self.attributes.get(transcript_id_key,'')
            if self.field[6] == '+':
                self.strand = 1
            elif self.field[6] == '-':
                self.strand = -1
            else:
                self.strand = 0
    
    def __eq__(self, other): return self.span == other.span
    def __ne__(self, other): return self.span != other.span
    def __gt__(self, other): return self.span >  other.span
    def __ge__(self, other): return self.span >= other.span
    def __lt__(self, other): return self.span <  other.span
    def __le__(self, other): return self.span <= other.span

cpdef RNAseqMapping anno_to_mapping_object(AnnotationObject parent, list children):
    """Given a top-level 'transcript' GTF/GFF feature and a list of 'exon' children,
    return a matching RNAseqMapping object."""
    pass

cpdef RNAseqMapping generate_read_from_gtf(RNAseqDataset dataset, list gtf_lines):
    """Converts a list of GTF-formatted strings to an RNAseqMapping object
    for a single transcript, retaining exon positional information and
    attributes of the 'transcript' parent line. Does not assume the lines are in
    any particular order, but does require that they belong to the same transcript.
    """
    pass

cpdef RNAseqMapping generate_read_from_gff(RNAseqDataset dataset, list gff_lines):
    """Converts a list of GFF3-formatted strings to an RNAseqMapping object
    for a single transcript, retaining exon positional information and
    attributes of the 'transcript' parent line. Does not assume the lines are in
    any particular order, but does require that they belong to the same transcript.
    """
    pass


##############################################################################################################
# BAM file processing #
##############################################################################################################
# A set of functions for generating an RNAseqMapping object from one or more lines of a SAM file
cpdef (int, int) range_of_reads(list reads):
    """Returns the leftmost and rightmost position of a 
    sorted list of RNAseqMapping objects."""
    cdef:
        RNAseqMapping read
        int leftmost, rightmost
    
    leftmost, rightmost = reads[0].span
    for read in reads:
        if read.span[1] > rightmost:
            rightmost = read.span[1]
    
    return (leftmost, rightmost)

cpdef np.ndarray calculate_coverage(list reads, int left, int right):
    """Returns a numpy array of coverage depth for a list of reads."""
    cdef:
        RNAseqMapping read
        int l, r
        (int, int) span
        np.ndarray coverage
        Py_ssize_t array_length

    # Iterate over reads to populate a coverage numpy array
    array_length = right-left
    coverage = np.zeros(array_length, dtype=np.float32)
    for read in reads:
        for span in read.ranges:
            l = span[0]-left
            r = span[1]-left
            coverage[l:r] += read.weight
    
    return coverage

def split_chunk(list reads, float minimum_proportion, int max_gap):
    """Takes a list of reads and cuts it at the first span of
    > max_gap nucleotides at a coverage < minimum_proportion * max_cov
    of the piece being traversed.
    """
    cdef:
        RNAseqMapping read
        int leftmost, rightmost, gap, number_of_reads, index, last_index
        (int, int) span
        np.ndarray coverage, is_gap
        Py_ssize_t array_length, i, l, r
        float c, max_cov
        bint was_gap
    
    # Iterate over reads to populate a coverage numpy array
    leftmost, rightmost = range_of_reads(reads)
    coverage = calculate_coverage(reads, leftmost, rightmost)
    array_length = rightmost-leftmost
    is_gap = np.zeros(array_length, dtype=np.bool)
    max_cov = 0
    # Iterate over the array, keeping track of places where a gap was found
    cdef float32[:] COV = coverage
    for i in range(array_length):
        c = COV[i]
        if c > max_cov:
            max_cov = c
            gap = 0
        elif c <= max_cov * minimum_proportion:
            gap += 1
            if gap > max_gap:
                if gap == max_gap + 1: # Gap just began, backfill
                    is_gap[i-max_gap:i] = True
                
                is_gap[i] = True
        else:
            if gap > 0: # A gap just terminated
                max_cov = c
            
            gap = 0
    
    # A final iteration over reads gives our generator:
    # IF the read is on the leading edge (extends rightmost)
    # AND IF the read is unspliced (no True in read.splice)
    # AND IF the read encounters a gap
    # THEN yield the set of reads up to (and not including)
    number_of_reads = len(reads)
    leading_edge = l = r = -1
    last_index = 0
    was_gap = False
    for index in range(number_of_reads): 
        read = reads[index]
        if read.span[1] > leading_edge: # Extends leading edge
            leading_edge = read.span[1]
            if not True in read.splice: # No splice junctions
                span = read.ranges[-1]
                l = span[0] - leftmost
                r = span[1] - leftmost
                if np.all(is_gap[l:r]): # The entire read is in a gap
                    if not was_gap: # First time encountering the gap
                        yield reads[last_index:index]
                    
                    last_index = index
                    was_gap = True
                else:
                    was_gap = False
    
    yield reads[last_index:]

cdef str get_flank(dict genome, str chrom, int pos, int strand, str label_type, int label_len):
    """Gets flanking region to a genomic position to compare to an end label
    """
    cdef str flank
    cdef int range_start, range_end
    flank = 'X'*label_len
    if (label_type == 'S' and strand == 1) or (label_type == 'E' and strand == -1):
        range_start = pos - label_len
        range_end = pos
    else:
        range_start = pos + 1
        range_end = pos + 1 + label_len
    
    flank = genome[chrom][range_start:range_end]
    if strand == -1:
        flank = fu.rc(flank)
    
    return flank

cdef int get_junction_strand(dict genome, str chrom, int left, int right):
    """ Returns 1(+), -1(-), or 0(.) for a left/right
    pair of splice junction positions based
    on the flanking genomic sequence """
    cdef str flanking_sequence
    cdef int strand
    flanking_sequence = get_flank(genome, chrom, left-1, 1, 'E', 2) + get_flank(genome, chrom, right, 1, 'S', 2)
    flanking_sequence = flanking_sequence.upper()
    if flanking_sequence in ['GTAG','GCAG','ATAC','GTGG']:
        strand = 1
    elif flanking_sequence in ['CTAC','CTGC','GTAT','CCAC']:
        strand = -1
    else:
        strand = 0
    
    return strand

# cdef (int, int) shift_junction(dict genome, str chrom, int left, int right, int junction_strand, int sj_shift):
#     """Given a (left, right) pair of splice junction positions,
#     searches for a canonical splice junction at the position. If one is not found,
#     performs a search in range +-sj_shift around the given positions for a canonical junction."""
#     cdef:
#         int j0, j1
#         str leftseq, rightseq
    
#     j0 = left
#     j1 = right
#     if strand == 1:
#         leftseq = 'GT'
#         rightseq = 'AG'
#     else:
#         leftseq = 'CT'
#         rightseq = 'AC'
    
#     left_flank = get_flank(genome, chrom, j0-1, 1, 'E', 2).upper()
#     right_flank = get_flank(genome, chrom, j1, 1, 'S', 2).upper()
#     if left_flank != leftseq:
#         strand = 1
#     elif flanking_sequence in ['CTAC','CTGC','GTAT','CCAC']:
#         strand = -1
#     else:
#         strand = 0
    
#     return (j0, j1)


def parse_MD_string(str mdstring):
    """Creates a generator object to yield the elements
    in an MD:Z tag from the SAM format. The elements are either
    one or more numbers or one or more characters."""
    cdef:
        bint was_digit, is_digit
        str current, c
    
    was_digit = True
    current = ''
    for c in mdstring:
        is_digit = c.isdigit()
        if is_digit != was_digit: # A transition happened
            yield current
            current = ''
            was_digit = is_digit
        
        current += c
    
    yield current



cpdef parse_SAM_CIGAR(int pos, list cigartuples, str mdstring, float error_rate=0.1):
    """Converts the po and pysam-parsed CIGAR object of a SAM file to 
    an array of 0-indexed open (left,right) blocks as ranges and gaps.
    CIGAR operators:
    M  0  Match; can be either an alignment match or mismatch
    I  1  Insertion; present in the read  but not in the reference
    D  2  Deletion; present in the reference but not in the read
    N  3  Skipped region; a region is not present in the read
    S  4  Soft Clipping;  the clipped nucleotides are present in the read
    H  5  Hard Clipping; not present in the read
    P  6  Padding; padded area in the read and not in the reference
    """
    cdef list ranges, gaps
    cdef int cigar_len, head, tail, current_pos, operator, o_len, i, leftside, rightside, match_count
    cdef bint first_element, jumped
    cdef str next_match
    ranges = []
    gaps = []
    head = 0
    tail = 0
    current_pos = pos
    first_element = True
    jumped = True
    cigar_len = len(cigartuples)
    mismatches = 0
    exon_size = 0
    i = 0
    mismatch_generator = parse_MD_string(mdstring)
    match_count = int(next(mismatch_generator))
    for i in range(cigar_len):
        operator = cigartuples[i][0]
        o_len = cigartuples[i][1]
        if operator == 4 or operator == 5: # Softclipped ranges are stored in either head or tail
            if first_element:
                head = o_len
            else:
                tail = o_len
        
        if operator == 0: # Match
            if jumped: # The last match was across an intron, dump it
                jumped = False
                if exon_size > 0 and mismatches/exon_size > error_rate: # Last exon didn't pass quality threshold
                    del ranges[-1]
                    if len(ranges) == 0:
                        head = -1
                
                ranges += [(current_pos,current_pos+o_len)]
                mismatches = 0
                exon_size = 0
            else: # Continuation of exon, update right side
                ranges[-1] = (ranges[-1][0], current_pos+o_len)
            
            exon_size += o_len
            match_count -= o_len
            current_pos += o_len
        elif operator == 3: # Skipped region (N)
            leftside = current_pos
            current_pos += o_len
            rightside = current_pos
            gaps += [(leftside,rightside)]
            jumped = True
        elif operator == 2: # Deletion in query
            current_pos += o_len
            mismatches += o_len
        elif operator == 1: # Insertion in query
            mismatches += o_len
        
        while match_count < 0:
            try:
                next_match = next(mismatch_generator)
                if next_match.isdigit():
                    match_count += int(next_match)
                elif next_match[0] != '^':
                    match_count += 1
                    mismatches += 1
            except StopIteration:
                print("# RAN OUT")
                match_count = 0

        if first_element:
            first_element = False
    
    if exon_size > 0 and mismatches/exon_size > error_rate: # Last exon didn't pass quality threshold
        del ranges[-1]
        tail = -1
    
    return ranges, gaps, head, tail

cdef (bint, bint, int, int) parse_tag(str string, str tagsplit='_TAG='):
    """Updates readtype based on a tag present in the ID string"""
    cdef int e_len, s_len
    cdef bint s_tag, e_tag
    cdef str ID_string
    cdef list ID_length
    s_len = 0
    e_len = 0
    s_tag = False
    e_tag = False
    if tagsplit in string: # If _TAG= exists in the read name, add these TAGS to the readtype
        ID_string = string.split(tagsplit)[-1].upper()
        if ID_string != '':
            ID_length = [int(i) for i in ID_string[1:].split('E') if i.isdigit()]
    else:
        ID_string = ''
        ID_length = []
    
    if ID_string != '':
        if 'S' in ID_string or 'C' in ID_string:
            s_tag = True
            s_len = ID_length[0]
            if 'E' in ID_string:
                e_tag = True
                e_len = ID_length[-1]
        elif 'E' in ID_string:
            e_tag = True
            e_len = ID_length[0]
    
    return s_tag, e_tag, s_len, e_len

cdef list generate_read_from_bam(RNAseqDataset dataset, list input_lines, bint ignore_ends=False, bint secondary=False):
    """Convert a list of pysam.AlignedSegment objects into RNAseqMappings that can be added to an RNAseqDataset"""
    cdef:
        int s_len, e_len, Nmap, counter, input_len, map_number, mate, strand
        int i, gap_len, pos, junction_strand, chrom_id, start_pos, end_pos, trim_pos
        str ID, chrom, js, seq, trimmed_nuc
        (bint, bint, int, int) ID_tags = (False, False, 0, 0)
        dict mappings
        list splice, gaps, ranges, introns
        bint is_upstream_side, is_downstream_side, stranded, stranded_method
        (int, int) g
        array.array flankmatch
    
    stranded = stranded_method = False
    if dataset.stranded: # The read is strand-specific
        stranded_method = True # The method to generate the read is inherently stranded
    
    input_len = len(input_lines)
    mappings = {} # Make an empty array of length to store each mapping
    s_len = len(dataset.start_array)
    e_len = len(dataset.end_array)
    Nmap = 0
    map_number = 0
    counter = 1
    ID = 'none'
    seq = ''
    cdef float weight = float(1)/input_len
    capped = dataset.capped
    for line in input_lines: # Each line must be a pysam.AlignedSegment
        if ID == 'none':
            ID_tags = (False, False, 0, 0)
            ID = line.query_name
            if not ignore_ends:
                ID_tags = parse_tag(ID)
            
            if ID_tags[2] > 0:
                if ID_tags[2] < s_len:
                    s_len = ID_tags[2]
            
            if ID_tags[3] > 0:
                if ID_tags[3] < e_len:
                    e_len = ID_tags[3]
        
        s_tag = dataset.s_tag or ID_tags[0]
        e_tag = dataset.e_tag or ID_tags[1]
        if s_tag or e_tag:
            stranded = True
        
        if Nmap == 0: # Update the mapping number with attribute NH:i
            if line.has_tag('NH'):
                Nmap = line.get_tag('NH')
        
        if line.has_tag('HI'): # Get which of Nmap mappings this line belongs to
            map_number = line.get_tag('HI') - 1
        else:
            map_number = counter - 1 
        
        if line.is_unmapped or line.is_supplementary: # Skip unmapped reads and poorly mapped reads
            counter += 1
            continue
        
        if line.is_secondary and not secondary: # Unless secondary alignments are allowed, skip these too
            counter += 1
            continue

        if seq == '':
            seq = line.query_sequence
        
        # Determine the RNA strand
        mate = 1
        strand = 0
        if line.is_paired:
            if not line.is_proper_pair:
                continue
            
            if line.is_read1:
                mate = 1
                if Nmap < counter:
                    Nmap = counter
                
                counter += 1
                if stranded:
                    if line.is_reverse:
                        strand = -1
                    else:
                        strand = 1
            else:
                mate = 2
                if stranded:
                    if line.is_reverse:
                        strand = 1
                    else:
                        strand = -1
        else:
            mate = 1
            if Nmap < counter:
                Nmap = counter
            
            counter += 1
            if stranded:
                if line.is_reverse:
                    strand = -1
                else:
                    strand = 1

        pos = line.reference_start
        chrom_id = line.reference_id
        try:
            chrom = line.header.get_reference_name(chrom_id)
        except:
            chrom = line.reference_name
        
        # Parse the SAM CIGAR string to get mapped positions, splice junction sites, and softclipped positions
        mdstring = line.get_tag('MD') if line.has_tag('MD') else ''
        ranges, introns, head, tail = parse_SAM_CIGAR(pos, line.cigartuples, mdstring, error_rate=dataset.mismatch_rate)
        if len(ranges) == 0: # No exons of passing quality were found
            continue

        # Check which gaps between exon blocks are present in intron blocks
        gaps = [(a,b) for a,b in zip([r for l,r in ranges[:-1]],[l for l,r in ranges[1:]])] # List of all gaps in ranges
        gap_len = len(gaps)
        splice = [False]*gap_len # List of booleans indicating whether each gap is a splice junction
        junction_strand = 0
        for i in range(gap_len):
            junction_strand = 0
            g = gaps[i]
            if g in introns:
                splice[i] = True
                if junction_strand == 0: # Try to resolve a nonstranded read with splice junctions
                    if line.has_tag('XS'):
                        js = line.get_tag('XS')
                        if js == '+':
                            junction_strand = 1
                        elif js == '-':
                            junction_strand = -1
                    else:
                        if dataset.genome:
                            junction_strand = get_junction_strand(dataset.genome, chrom, g[0], g[1])
                            # if junction_strand == 0: # Did not find a valid splice junction
                            #     j0, j1, s = shift_junction(genome, chrom, g[0], g[1], sj_shift)
                            #     if j0 == -1 or j1 == -1:
                            #         splice[i] = False
                            #     else:
                            #         # If shift_junction() found a shift, update the appropriate edge in ranges
                            #         junction_strand = s
                            #         if j0 != g[0]:
                            #             ranges[i][1] = j0
                            #         if j1 != g[1]:
                            #             ranges[i+1][0] = j1
            
                if junction_strand == 0:
                    splice[i] = False
        
        # Quality control of end labels:
        # 1) An improperly mapped 5'/3' end of a read should be stripped of its tag
        # 2) Upstream untemplated Gs on a 5' end are evidence for a cap structure (requires genome)
        # 3) End labels that match the genome-templated sequence are false positive trims (requires genome)
        # 4) False positive oligo-dT priming can occur at genome-templated purine-rich sites (requires genome)
        # 5) False positive template-switching can occur at matching RNA sites of 3+ nucleotides (requires genome)
        
        # EVALUATE SOFTCLIPPED NUCLEOTIDES
        is_upstream_side = mate == 1
        is_downstream_side = (mate == 1 and not line.is_paired) or mate == 2
        if head == -1: # Special case: left side clipped off for quality issues
            if strand == 1:
                s_tag = capped = False
            elif strand == -1:
                e_tag = False
        
        if tail == -1: # Special case: right side clipped off for quality issues
            if strand == 1:
                e_tag = False
            elif strand == -1:
                s_tag = capped = False
        
        if s_tag: # Check for upstream untemplated Gs
            if is_upstream_side:
                if strand == 1:
                    if head > 0: # Plus-stranded left clip
                        if head > 4: # Softclipped sequence is too long
                            s_tag = False
                            capped = False
                        else: # From 1-4 softclipped nucleotides; check if untemplated G's
                            if seq[:head] == 'G'*head: # Softclipped nucleotides are G
                                if dataset.genome:
                                    if get_flank(dataset.genome, chrom, ranges[0][0], 1, 'S', 1) != 'G': # The flanking nucleotide is NOT G
                                        capped = True # One or more upstream untemplated Gs were detected
                                else:
                                    capped = True
                elif strand == -1:
                    if tail > 0: # Minus-stranded right clip
                        if tail > 4: # Softclipped sequence is too long
                            s_tag = False
                            capped = False
                        else:
                            if seq[-tail:] == 'C'*tail: # Sofclipped nucleotides are (antisense) G
                                if dataset.genome:
                                    if get_flank(dataset.genome, chrom, ranges[-1][-1]-1, -1, 'S', 1) != 'G': # The flanking nucleotide is NOT G
                                        capped = True
                                else:
                                    capped = True
            elif is_downstream_side: # Can't have an s_tag if downstream mate
                s_tag = False
                capped = False
                
        
        if e_tag: # Check for softclipped nucleotides to add back
            if is_downstream_side:
                if strand == 1:
                    if tail > 0:
                        if tail > 10: # Softclipped sequence is too long
                            e_tag = False
                        elif tail < 4:
                            ranges[-1] = (ranges[-1][0],ranges[-1][1]+tail)
                elif strand == -1:
                    if head > 0: # Left clip exists
                        if head > 10: # Left clip is too long
                            e_tag = False
                        elif head < 4:
                            ranges[0] = (ranges[0][0]-head,ranges[0][1])
            else: # Can't have an s_tag if it is the upstream mate
                e_tag = False
        
        # EVALUATE FALSE POSITIVE TAGS
        start_pos = 0
        end_pos = 0
        if strand == 1:
            start_pos = ranges[0][0]
            end_pos = ranges[-1][-1]-1
        elif strand == -1:
            start_pos = ranges[-1][-1]-1
            end_pos = ranges[0][0]
        
        if dataset.genome:            
            if s_tag and is_upstream_side:
                flank = get_flank(dataset.genome, chrom, start_pos, strand, 'S', s_len) # Get upstream flanking sequence to start
                if len(flank) > 0:
                    flankmatch = dataset.start_array[-s_len:]
                    if fu.IUPACham(fu.nuc_to_int(flank), flankmatch, dataset.mismatch_rate*s_len) <= dataset.mismatch_rate*s_len:
                        s_tag = False # Query sequence matched well enough to masking sequence
            
            if e_tag and is_downstream_side:
                flank = get_flank(dataset.genome, chrom, end_pos, strand, 'E', e_len) # Get downstream flanking sequence to end
                if len(flank) > 0:
                    flankmatch = dataset.end_array[:e_len]
                    if fu.IUPACham(fu.nuc_to_int(flank), flankmatch, dataset.mismatch_rate*e_len) <= dataset.mismatch_rate*e_len:
                        e_tag = False # Query sequence matched well enough to masking sequence
        
        # Reconcile strand information given by start, end, and splice
        if junction_strand != 0:
            if strand != junction_strand: # Splice disagrees with end tags; remove tags
                strand = junction_strand
                s_tag = False
                e_tag = False
                capped = False
        
        if not stranded_method and not s_tag and not e_tag and junction_strand == 0:
            strand = 0 # No strand information can be found
        
        # Generate a ReadObject with the parsed attributes above
        read_data = ELdata(chrom_id, 0, strand, ranges, splice, s_tag, e_tag, capped, round(weight,2))
        current_mapping = RNAseqMapping(read_data)
        
        if map_number not in mappings:
            mappings[map_number] = current_mapping
        else:
            mappings[map_number].merge(current_mapping) # merge two mate-pair ReadObjects together
    
    return list(mappings.values())


def read_generator(fileconn, RNAseqDataset dataset, str file_type, int max_gap):
    """Yields a contiguous chunk of reads from the input file
    separated on either side by a gaps > max_gap"""
    cdef RNAseqMapping last_read
    cdef int rightmost, last_chrom
    if file_type == 'elr':
        add_read = dataset.add_read_from_ELR
    elif file_type == 'bed':
        add_read = dataset.add_read_from_BED
    elif file_type == 'bam' or file_type == 'sam':
        add_read = dataset.add_read_from_BAM
    else:
        return
    
    rightmost = -1
    last_chrom = -1
    for line in fileconn:
        if type(line) is str:
            if line[0] == '#':
                header_line = line.rstrip().split(' ')
                if header_line[0] == '#S':
                    dataset.add_source(header_line[-1])
                elif header_line[0] == '#C':
                    dataset.add_chrom(header_line[-1])

                continue
        
        add_read(line)
        last_read = dataset.read_list[-1]
        if rightmost > 0:
            if last_read.left() > (max_gap + rightmost): # A sufficiently large gap was jumped
                yield dataset.read_list[:-1]
                dataset.read_list = [last_read]
            elif last_read.chrom != last_chrom:
                yield dataset.read_list[:-1]
                dataset.read_list = [last_read]
                rightmost = -1
        
        if last_read.span[1] > rightmost:
            rightmost = last_read.right()
        
        last_chrom = last_read.chrom

    yield dataset.read_list
    fileconn.close()


##############################################################################################################
##############################################################################################################
# Object for summarizing the coverage and features of a collection of RNAseqMapping objects #
##############################################################################################################
##############################################################################################################

cdef class BranchpointArray:
    cdef readonly dict S_plus, S_minus, E_plus, E_minus, J_plus, J_minus
    cdef readonly int leftmost, rightmost, extend, end_extend, length, min_overhang
    cdef readonly tuple branchpoints
    cdef readonly float weight, threshold, cap_percent, minimum_proportion
    cdef readonly np.ndarray depth
    cdef public OrderedArray bp_plus, bp_minus
    def __init__(self, int leftmost, int rightmost, tuple reads, int extend, int end_extend, float minimum_proportion, float cap_percent=0.0, int min_overhang=0):
        """Makes a collection of positions in the locus that act as dividing
        points for all nodes of the graph. Every donor (D) and acceptor (A) 
        site are retained as-is, but start (S) and end (E) are collapsed into
        a set of reference points."""
        cdef BranchPoint BP, bp, lbp, rbp, lterm, rterm
        cdef RNAseqMapping read
        cdef float e_threshold, s_threshold, s_counter, e_counter, weight, wt, threshold_depth
        cdef str branchtype, junction_hash
        cdef int pos, dist, bt, length, i, number_of_reads
        cdef char strand
        cdef list s_rank, e_rank, rank_sort, merged_plus, merged_minus, to_remove, donors, acceptors, gaps, gap_branchpoints
        cdef (float, int, int, int) ranking
        cdef (int, float) item
        cdef (int, int) block
        cdef dict lookup, end_counts
        cdef set donor_sites, acceptor_sites
        cdef bint first_element, s_added
        
        Sp, Sm, Cp, Cm, Ep, Em, Dp, Dm, Ap, Am, S, E, bp_dict = Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), Counter()
        self.J_plus = {}
        self.J_minus = {}
        self.extend = extend
        self.end_extend = end_extend
        self.weight = 0
        self.minimum_proportion = minimum_proportion
        self.min_overhang = min_overhang
        self.threshold = 1 - self.minimum_proportion
        self.cap_percent = cap_percent
        self.leftmost = leftmost
        self.rightmost = rightmost
        length = self.rightmost-self.leftmost
        self.depth = calculate_coverage(list(reads), self.leftmost, self.rightmost)
        self.weight = np.sum(self.depth)
        number_of_reads = len(reads)
        for i in range(number_of_reads):
            read = reads[i]
            if read.strand == 1:
                for l,r in read.junctions():
                    Dp[l] += read.weight
                    Ap[r] += read.weight
                    block = (l,r)
                    junction_hash = str(block)
                    self.J_plus[junction_hash] = self.J_plus.get(junction_hash, 0) + read.weight
                
                if read.s_tag:
                    pos = read.span[0]
                    Sp[pos] += read.weight
                    if read.capped:
                        Cp[pos] += read.weight
                
                if read.e_tag:
                    pos = read.span[1]
                    Ep[pos] += read.weight
            elif read.strand == -1:
                for l,r in read.junctions():
                    Am[l] += read.weight
                    Dm[r] += read.weight
                    block = (l,r)
                    junction_hash = str(block)
                    self.J_minus[junction_hash] = self.J_minus.get(junction_hash, 0) + read.weight
                
                if read.e_tag:
                    pos = read.span[0]
                    Em[pos] += read.weight
                
                if read.s_tag:
                    pos = read.span[1]
                    Sm[pos] += read.weight
                    if read.capped:
                        Cm[pos] += read.weight
        
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
                (-weight, self.distance_from_edge(pos, strand, 'S'), 2, pos)
                for pos, weight in S.most_common()
            ]
            e_rank = [
                (-weight, self.distance_from_edge(pos, strand, 'E'), 1, pos)
                for pos, weight in E.most_common()
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
                        BP = StartPoint(strand, pos, weight, self.end_extend, C[pos])
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

        merged_plus = []
        first_element = True
        for bp in self.bp_plus:
            if bp.weight >= self.minimum_proportion * end_counts[bp.branchtype]:
                if bp.branchtype != 'S': # No extra restrictions for non-S branchpoints
                    merged_plus.append(bp)
                elif first_element:  # Is the leftmost branchpoint
                    merged_plus.append(bp)
                elif bp.percent_capped() >= self.cap_percent: # Passes the cap threshold
                    merged_plus.append(bp)
                
                first_element = False       

        end_counts = {branchtype:0.0 for branchtype in ['S','E','A','D','N']}
        for bp in self.bp_minus:
            if bp.branchtype == 'S':
                end_counts['S'] += bp.weight
            elif bp.branchtype == 'E':
                end_counts['E'] += bp.weight

        merged_minus = []
        for  bp in self.bp_minus:
            if bp.weight >= self.minimum_proportion * end_counts[bp.branchtype]:
                if bp.branchtype != 'S': # No extra restrictions for non-S branchpoints
                    merged_minus.append(bp)
                elif bp.percent_capped() >= self.cap_percent: # Passes the cap threshold
                    merged_minus.append(bp)
        
        # Allow the rightmost S- to be added regardless of cap_percent
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


        self.branchpoints = tuple(sorted(merged_plus + merged_minus))
        # Evaluate whether terminal branchpoints must be added.
        # If no boundary Start/Endpoint exists, add a nonspecified one at leftmost/rightmost
        passes_threshold = np.where(self.depth >= threshold_depth)[0]
        if passes_threshold.shape[0] > 0:
            lbp = BranchPoint('N', 0, self.leftmost+passes_threshold[0], 0)
            rbp = BranchPoint('N', 0, self.leftmost+passes_threshold[-1], 0)
        else: # No positions pass the threshold
            lbp = BranchPoint('N', 0, self.leftmost, 0)
            rbp = BranchPoint('N', 0, self.rightmost, 0)
        
        if len(self.branchpoints) == 0:
            self.branchpoints = tuple([lbp, rbp])
        else:
            # Extend left border if S>/E< doesn't contain
            lterm = self.branchpoints[0]
            rterm = self.branchpoints[-1]
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
                    self.branchpoints = tuple([lbp] + list(self.branchpoints) + [rbp])
                elif add_lbp:
                    self.branchpoints = tuple([lbp] + list(self.branchpoints))
                else:
                    self.branchpoints = tuple(list(self.branchpoints) + [rbp])
    
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
    
    cpdef list bp_from_counter(self, counter, str branchtype, char strand):
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


bp_typeorder = {'N':-1, 'E':0, 'A':1, 'D':2, 'S':3} # Sort order for branchpoint types
cdef class BranchPoint:
    """Represents a reference point for a Start, End, Donor, or Acceptor site."""
    cdef readonly str branchtype
    cdef readonly char strand
    cdef public int pos, index, left, right
    cdef public float weight
    cdef public (int, int) comparator
    def __init__(self, str branchtype, char strand, int pos, float weight):
        assert strand >= -1 and strand <= 1
        self.strand = strand
        self.branchtype = branchtype
        self.pos = pos
        self.left = self.right = self.pos
        self.weight = float(weight)
        if self.strand == 1:
            self.comparator = (self.pos, bp_typeorder[self.branchtype])
        else:
            self.comparator = (self.pos, -bp_typeorder[self.branchtype])
    
    cpdef float percent_capped(self):
        return 1.0

    def __repr__(self):
        strand = ['.','+','-'][self.strand]
        return '{}{}{} ({})'.format(self.branchtype, strand, self.pos, self.weight)
    
    def __eq__(self, other): return self.comparator == other.comparator
    def __ne__(self, other): return self.comparator != other.comparator
    def __gt__(self, other): return self.comparator >  other.comparator
    def __ge__(self, other): return self.comparator >= other.comparator
    def __lt__(self, other): return self.comparator <  other.comparator
    def __le__(self, other): return self.comparator <= other.comparator
    
    def write_as_bed(self, chrom):
        print('{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, self.left, self.right, self.index, self.weight, {-1:'-',1:'+',0:'.'}[self.strand]))


cdef class StartPoint(BranchPoint):
    """Represents a transcript start site"""
    cdef public int extend
    cdef public object span
    cdef public float capped
    def __init__(self, strand, pos, weight, extend=0, capped_weight=0):
        BranchPoint.__init__(self, 'S', strand, pos, weight)
        self.extend = extend
        self.capped = capped_weight
        self.span = range(self.left-self.extend, self.right+self.extend+1)
    
    cpdef float percent_capped(self):
        return self.capped/self.weight
    
    def merge(self, StartPoint other):
        """Merge other into self"""
        self.weight += other.weight
        self.capped += other.capped
        self.left = min(self.left, other.left)
        self.right = max(self.right, other.right)
        self.span = range(self.left-self.extend, self.right+self.extend+1)

cdef class EndPoint(BranchPoint):
    """Represents a transcript termination site"""
    cdef public int extend
    cdef public object span
    def __init__(self, strand, pos, weight, extend=0):
        BranchPoint.__init__(self, 'E', strand, pos, weight)
        self.extend = extend
        self.span = range(self.left-self.extend, self.right+self.extend+1)

    def merge(self, EndPoint other):
        """Merge other into self"""
        self.weight += other.weight
        self.left = min(self.left, other.left)
        self.right = max(self.right, other.right)
        self.span = range(self.left-self.extend, self.right+self.extend+1)

cdef class DonorPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight):
        BranchPoint.__init__(self, 'D', strand, pos, weight)

cdef class AcceptorPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight):
        BranchPoint.__init__(self, 'A', strand, pos, weight)

cdef class OrderedArray:
    """Maintains an array in sort order and a lookup dict for the location of each
    item that was added."""
    cdef readonly int n
    cdef readonly list ordered_array
    def __init__(self, array=[]):
        self.n = len(array)
        self.ordered_array = sorted([(v,i) for i,v in enumerate(array)])
    
    def __repr__(self):
        return 'OA({})'.format(self.n)
    
    def ordered_index(self, position=None):
        if position is None:
            return [i for v,i in self.ordered_array]
        else:
            return self.ordered_array[position][-1]
    
    def ordered_values(self, position = None):
        if position is None:
            return [v for v,i in self.ordered_array]
        else:
            return self.ordered_array[position][0]
    
    def __iter__(self):
        return iter(self.ordered_values())
    
    def __len__(self):
        return len(self.ordered_array)
    
    def __add__(self,other):
        new_OA = copy.deepcopy(self)
        for item in other:
            new_OA.add(item)
        
        return new_OA
    
    cdef int get_insert_pos(self, search_value, init=None):
        """Runs a binary search, returns an insert position where 'search_value'
        fits in array. Assumes array to be sorted. In the case of a tie,
        priority goes to the existing item."""
        cdef int insert_pos, first, last, mid
        insert_pos = -1
        if init is None:
            first = 0
            last = len(self.ordered_array) - 1
        else: # Start at an initialized index position (if a good guess of the correct index exists)
            if search_value == self.ordered_array[init]: # Exact match, insert pos is right of init
                return init + 1
            elif search_value > self.ordered_array[init]: # insert pos must be right of init
                first = init + 1
                last = len(self.ordered_array) - 1
            else: # insert pos must be left of init
                first = 0
                last = init - 1
        
        while first <= last and insert_pos == -1:
            mid = (first + last)//2 # Get the midpoint of the search space
            if search_value == self.ordered_array[mid]: # Exact match, insert pos is right of mid
                insert_pos = mid + 1
            elif search_value > self.ordered_array[mid]: # insert pos must be right of the midpoint
                first = mid + 1
            else: # insert pos must be left of the midpoint
                last = mid - 1
        
        if insert_pos == -1: # An exact match wasn't found, but the search converged
            insert_pos = first
        
        return insert_pos
    
    def add(self, item, init=None):
        """Add an item to the ordered array."""
        self.n += 1
        ipos = self.get_insert_pos((item, self.n), init)
        self.ordered_array.insert(ipos, (item, self.n))
    
    def delete(self, index):
        """Deletes the item at index position in the ordered array."""
        self.n -= 1
        del self.ordered_array[index]
    
    def add_list(self, list_to_add):
        for item in list_to_add:
            self.add(item)
    
    def probe(self, item, init=None):
        """Return the insert position of an item if it were to be added to the array"""
        return self.get_insert_pos((item, self.n + 1), init)

cpdef list get_gaps(np.ndarray[float, ndim=1] array, int maxgap, threshold = float(1)):
    """Returns True if no gap longer than maxgap in the array falls below the threshold"""
    cdef Py_ssize_t i, l
    cdef int gap_length, gap_left, gap_right
    cdef float value
    cdef bint gap_is_long_enough = False
    cdef list gaps = []
    cdef float [:] ARRAY = array
    gap_length = 0
    l = ARRAY.shape[0]
    for i in range(l):
        value = ARRAY[i]
        if value < threshold: # A gap begins or extends
            if gap_length == 0: # Beginning of a gap
                gap_left = i
            
            gap_length += 1
            if not gap_is_long_enough:
                if gap_length > maxgap:
                    gap_is_long_enough = True
        else: # Not in a gap
            if gap_is_long_enough: # The last position was the end of a gap
                gap_right = i
                gaps.append((gap_left, gap_right))
            
            gap_length = 0
            gap_is_long_enough = False
    
    return gaps
