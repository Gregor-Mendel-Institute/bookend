#cython: language_level=3
cimport cython
import array
from cpython cimport array
import numpy as np
cimport numpy as np
import bookend.core.cython_utils._fasta_utils as fu
from collections import namedtuple, Counter
from ast import literal_eval
import copy
ctypedef unsigned char uint8
ctypedef np.float32_t float32

##############################################
# Object for defining an RNA sequencing read #
##############################################
ELdata = namedtuple('ELdata', 'chrom source strand ranges splice s_tag e_tag capped weight')
cdef class RNAseqMapping():
    cdef public int chrom, source, strand, s_len, e_len
    cdef public list ranges, splice
    cdef public bint s_tag, e_tag, capped, complete, is_reference
    cdef public dict attributes
    cdef public (int, int) span
    cdef public float weight, coverage
    def __init__(self, input_data, attributes = None):
        """Initializes a Read Object given a tuple of input data.
        Requires a chromosome, strand, source, weight, a sorted tuple of
        exon ranges and an array of booleans indicating which gaps between exons are splice junctions."""
        self.is_reference = False
        self.s_len = self.e_len = 0
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

    cpdef gaps(self):
        """Returns an array of 0-indexed (start, end) tuples of gaps between ranges"""
        if len(self.ranges) == 1:
            return []
        
        return [(self.ranges[i][-1], self.ranges[i+1][0]) for i in range(len(self.ranges)-1)]
    
    cpdef junctions(self):
        """Returns an array of 0-indexed (start, end) tuples of intron locations"""
        j_array = []
        for i,j in enumerate(self.splice):
            if j:
                j_array += [(self.ranges[i][-1], self.ranges[i+1][0])]
        
        return j_array
    
    cpdef int diff(self, other):
        """Returns the total number of nucleotides overlapped by only one read."""
        return abs(self.span[0]-other.span[0])+abs(self.span[1]-other.span[1])
    
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
    
    cpdef bint splice_match(self, RNAseqMapping other, bint ignore_ends=True):
        """Returns bool of whether self and other ranges match perfectly.
        Discards 5' and 3' terminus if ignore_ends."""
        cdef list exons_self, exons_other
        if ignore_ends:
            return self.junctions() == other.junctions()
        else:
            return self.ranges == other.ranges

    cpdef bint antisense_match(self, RNAseqMapping other, int end_extend):
        """Returns bool if self is antisense to other and overlaps by at least end_extend."""
        cdef int left, right, self_size, other_size
        if self.overlaps(other):
            self_size = self.span[1]-self.span[0]
            other_size = other.span[1]-other.span[0]
            end_extend = min([end_extend, self_size, other_size])
            if self.strand == -other.strand:
                left, right = self.overlap_range(other)
                if right - left >= end_extend:
                    return True
            
        return False
    
    cpdef bint sense_match(self, RNAseqMapping other, int end_extend):
        """Returns bool if self is sense to other and overlaps by at least end_extend."""
        cdef int left, right, self_size, other_size
        if self.overlaps(other):
            self_size = self.span[1]-self.span[0]
            other_size = other.span[1]-other.span[0]
            end_extend = min([end_extend, self_size, other_size])
            if self.strand == other.strand:
                left, right = self.overlap_range(other)
                if right - left >= end_extend:
                    return True
            
        return False

    cpdef bint is_compatible(self, RNAseqMapping other, bint ignore_ends=False, bint ignore_source=False):
        """Self and other contain no attributes that demonstrate they could
        not be subsequences of a common longer molecule."""
        cdef (int, int) overlap
        if self.chrom != other.chrom:
            return False
        
        if not ignore_source:
            if self.source != other.source:
                return False
        
        if self.strand != 0 and other.strand != 0 and self.strand != other.strand:
            return False
        
        # If self or other contain terminal tags, the tags must be compatible
        if not ignore_ends:
            if self.ends_clash(other):
                return False
        
        if not self.overlaps(other):
            return True # No incompatibilities were found
        
        # If the two reads share a chrom and strand and overlap,
        # check the overlapping range for identical splice architecture
        overlap = self.overlap_range(other)
        j1 = [j for j in self.junctions() if j[1] > overlap[0] and j[0] < overlap[1]]
        j2 = [j for j in other.junctions() if j[1] > overlap[0] and j[0] < overlap[1]]
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
    
    def get_node_labels(self, record_artifacts=False):
        """Returns a string with one label for each edge of each range in self.ranges."""
        if self.strand == 1:
            gapchar = 'DA'
            if self.s_tag:
                if self.capped:
                    startchar = 'C'
                else:
                    startchar = 'S'
            else:
                if record_artifacts and self.s_len > 0:
                    startchar = 's'
                else:
                    startchar = '.'
            
            if self.e_tag:
                endchar = 'E'
            else:
                if record_artifacts and self.e_len > 0:
                    endchar = 'e'
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
                if record_artifacts and self.s_len > 0:
                    endchar = 's'
                else:
                    endchar = '.'
            
            if self.e_tag:
                startchar = 'E'
            else:
                if record_artifacts and self.e_len > 0:
                    startchar = 'e'
                else:
                    startchar = '.'
        else:
            gapchar = '..'
            startchar = endchar = '.'
            if record_artifacts:
                if self.s_len > 0:
                    startchar = 's'
                
                if self.e_len > 0:
                    endchar = 'e'
        
        return(''.join([startchar]+[gapchar if i else '..' for i in self.splice]+[endchar]))
    
    cpdef write_as_elr(self, as_string=True, record_artifacts=False):
        """Returns a string that represents the ReadObject
        in the end-labeled read (ELR) format"""
        cdef str elr_strand, labels
        cdef list block_ends, elr_line
        elr_strand = '.'
        if self.strand == 1:
            elr_strand = '+'
        elif self.strand == -1:
            elr_strand = '-'
        
        block_ends = flatten(self.ranges)
        lengths = [block_ends[i]-block_ends[i-1] for i in range(1,len(block_ends))]
        labels = self.get_node_labels(record_artifacts)
        EL_CIGAR = ''.join([str(a)+str(b) for a,b in zip(labels,lengths+[''])])
        read_len = self.right() - self.left()
        elr_line = [self.chrom, self.left(), read_len, elr_strand, EL_CIGAR, self.source, round(self.weight,2)]
        if as_string:
            return '\t'.join([str(i) for i in elr_line])
        else:
            return elr_line
     
    cpdef write_as_bed(self, chrom_array, source_array, as_string=True, score_column='weight', record_artifacts=False, name_attr=None, color=None):
        """Returns a string that represents the ReadObject
        in a 15-column BED format"""
        labels = self.get_node_labels(record_artifacts)
        bed_strand = '.'
        if self.strand == 1:
            bed_strand = '+'
        elif self.strand == -1:
            bed_strand = '-'
        
        l = labels[0]
        r = labels[-1]
        ends = l+r
        if color is None:
            rgb = '0,0,0'
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
        else:
            rgb = color
        
        name = '.'
        if name_attr is not None:
            name = self.attributes.get(name_attr, '.')
        
        if score_column == 'weight':
            score = round(self.weight,2)
        elif score_column == 'coverage':
            score = round(self.coverage,2)
        else:
            score = str(score_column)
        
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


############################################################
# Object for defining a collection of RNA sequencing reads #
############################################################
config_defaults = {
    'source':'',
    's_tag':False,
    'e_tag':False,
    'capped':False,
    'stranded':False,
    'start_seq':'ACGGG',
    'end_seq':'RRRRRRRRRRRRRRR',
    'minlen_strict':20,
    'minlen_loose':25,
    'mismatch_rate':0.2,
    'min_reps':2,
    'cap_bonus':5,
    'confidence_threshold':0.5
}

cdef class RNAseqDataset():
    cdef public list read_list, chrom_array, source_array, chrom_lengths
    cdef public int chrom_index, source_index
    cdef public dict chrom_dict, source_dict
    cdef readonly dict config, genome, label_tally
    cdef readonly bint s_tag, e_tag, capped, stranded, ignore_ends
    cdef readonly str start_seq, end_seq
    cdef readonly int minlen, minlen_strict, minlen_loose
    cdef readonly float mismatch_rate
    cdef readonly array.array start_array, end_array

    def __init__(self, chrom_array=None, source_array=None, chrom_lengths=None, genome_fasta=None, config=config_defaults):
        """Container for RNAseqMapping objects. Stores a reference dictionary for all
        chromosome names and sample names. Contains methods for parsing
        a variety of files into a collection of read objects."""
        self.label_tally = {'S':Counter(), 's':Counter(), 'E':Counter(), 'e':Counter()}
        self.read_list = []
        self.config = config
        self.s_tag = self.config['s_tag']
        self.e_tag = self.config['e_tag']
        self.capped = self.config['capped']
        self.stranded = self.config['stranded']
        self.start_seq = self.config['start_seq']
        self.end_seq = self.config['end_seq']
        self.minlen_strict = self.config['minlen_strict']
        self.minlen_loose = self.config['minlen_loose']
        self.minlen = self.minlen_strict
        self.mismatch_rate = self.config['mismatch_rate']
        self.start_array = fu.nuc_to_int(self.start_seq)
        self.end_array = fu.nuc_to_int(self.end_seq)
        self.chrom_lengths = chrom_lengths
        self.chrom_dict = {}
        self.chrom_index = 0
        self.chrom_array = []
        if chrom_array is not None:
            for c in chrom_array:
                self.add_chrom(c)
        
        if genome_fasta is not None:
            self.genome, index = fu.import_genome(genome_fasta)
            if chrom_array is None:
                index_lines = [l.split('\t') for l in index.rstrip().split('\n')]
                self.chrom_array = [l[0] for l in index_lines]
                self.chrom_lengths = [int(l[1]) for l in index_lines]
                self.chrom_index = len(self.chrom_array)
                self.chrom_dict = dict(zip(self.chrom_array, range(self.chrom_index)))
        else:
            self.genome = {}
        
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
        cdef RNAseqMapping read
        cdef BAMobject BAM
        if type(bam_lines) is not list:
            bam_lines = [bam_lines]
        
        BAM = BAMobject(self, bam_lines, ignore_ends, secondary)
        new_read_list = BAM.generate_read()
        if len(new_read_list) > 0:
            read = new_read_list[0]
            if read.s_len > 0:
                if read.s_tag:
                    self.label_tally['S'][read.s_len] += 1
                else:
                    self.label_tally['s'][read.s_len] += 1
            
            if read.e_len > 0:
                if read.e_tag:
                    self.label_tally['E'][read.e_len] += 1
                else:
                    self.label_tally['e'][read.e_len] += 1
        
        self.read_list += new_read_list

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


######################################
# Utilities for GTF/GFF file parsing #
######################################
gtf_defaults = {
    'parent_types':set(['transcript']),
    'parent_key_transcript':'transcript_id',
    'parent_key_gene':'gene_id',
    'child_types':set(['exon']),
    'child_key_transcript':'transcript_id',
    'child_key_gene':'gene_id'
}
gff_defaults = {
    'parent_types':set([
        'mRNA','transcript',
        'snoRNA','tRNA','snRNA','rRNA','miRNA','ncRNA','mRNA_TE_gene','pseudogenic_transcript',
        'antisense_lncRNA','antisense_RNA','lnc_RNA']),
    'parent_key_transcript':'ID',
    'parent_key_gene':'Parent',
    'child_types':set(['exon','pseudogenic_exon']),
    'child_key_transcript':'Parent',
    'child_key_gene':'gene_id'
}
gtf_colorcode = {
    '3prime_overlapping_ncRNA': '203,98,130',
    'antisense': '194,71,71', 
    'antisense_RNA': '194,71,71', 
    'atlncRNA': '56,114,168', 
    'atRNA': '56,114,168', 
    'bidirectional_promoter_lncRNA': '153,108,206', 
    'IG_C_gene': '110,66,19', 
    'IG_C_pseudogene': '80,80,80', 
    'IG_D_gene': '110,66,19', 
    'IG_D_pseudogene': '80,80,80', 
    'IG_J_gene': '110,66,19', 
    'IG_LV_gene': '110,66,19', 
    'IG_pseudogene': '80,80,80', 
    'IG_V_gene': '110,66,19', 
    'IG_V_pseudogene': '80,80,80', 
    'lincRNA': '56,114,168', 
    'lncRNA': '56,114,168', 
    'macro_lncRNA': '6,74,140', 
    'miRNA': '198,95,84', 
    'misc_RNA': '249,185,54', 
    'Mt_rRNA': '0,0,0', 
    'Mt_tRNA': '0,0,0', 
    'ncRNA': '249,185,54',
    'nonsense_mediated_decay': '180,155,100', 
    'non_stop_decay': '180,155,100', 
    'nontranslating_CDS': '180,155,100',
    'otherRNA': '56,114,168',
    'polymorphic_pseudogene': '80,80,80', 
    'pre_miRNA': '198,95,84', 
    'processed_pseudogene': '80,80,80', 
    'processed_transcript': '180,155,100', 
    'protein_coding': '49,132,44', 
    'pseudogene': '80,80,80', 
    'retained_intron': '180,155,100', 
    'ribozyme': '249,185,54', 
    'rRNA': '0,0,0', 
    'RNase_MRP_RNA': '249,185,54', 
    'scaRNA': '249,185,54', 
    'scRNA': '249,185,54', 
    'sense_intronic': '56,114,168', 
    'sense_overlapping': '56,114,168', 
    'snoRNA': '249,185,54', 
    'snRNA': '249,185,54', 
    'SRP_RNA': '249,185,54', 
    'sRNA': '249,185,54', 
    'TEC': '120,120,120', 
    'transcribed_processed_pseudogene': '80,80,80', 
    'transcribed_unitary_pseudogene': '80,80,80', 
    'transcribed_unprocessed_pseudogene': '120,120,120', 
    'translated_processed_pseudogene': '80,80,80', 
    'translated_unprocessed_pseudogene': '120,120,120', 
    'tRNA': '0,0,0',
    'TR_C_gene': '110,66,19', 
    'TR_D_gene': '110,66,19', 
    'TR_J_gene': '110,66,19', 
    'TR_J_pseudogene': '80,80,80', 
    'TR_V_gene': '110,66,19', 
    'TR_V_pseudogene': '80,80,80', 
    'unitary_pseudogene': '80,80,80', 
    'unprocessed_pseudogene': '120,120,120'
}


cdef class AnnotationObject:
    cdef public dict attributes
    cdef public bint keep, parent
    cdef public str format, gene_id, transcript_id, chrom, source, anno_type
    cdef public int strand
    cdef public (int, int) span
    cdef public tuple fields
    def __init__(self, anno_string, format, config_dict):
        """Generates an intermediate object from a single line of a GTF/GFF file.
        This will not completely represent a transcript, only one exon."""
        cdef:
            set child_types, parent_types
            str gene_id_key, transcript_id_key
        
        anno_string = anno_string.rstrip()
        self.format = format
        self.keep = False
        self.parent = False
        self.gene_id = ''
        self.transcript_id = ''
        self.anno_type = ''
        self.fields = ()
        if len(anno_string) > 0:
            if anno_string[0] != '#':
                self.fields = tuple(anno_string.split('\t')[0:9])
                self.anno_type = self.fields[2]
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
                    self.attributes = self.parse_attributes(self.fields[8], self.format)
                    if self.parent:
                        gene_id_key=config_dict['parent_key_gene']
                        transcript_id_key=config_dict['parent_key_transcript']
                    else:
                        gene_id_key=config_dict['child_key_gene']
                        transcript_id_key=config_dict['child_key_transcript']
                    
                    self.gene_id = self.attributes.get(gene_id_key,'')
                    self.transcript_id = self.attributes.get(transcript_id_key,'')
                    self.attributes['gene_id'] = self.gene_id
                    self.attributes['transcript_id'] = self.transcript_id
                    if self.fields[6] == '+':
                        self.strand = 1
                    elif self.fields[6] == '-':
                        self.strand = -1
                    else:
                        self.strand = 0
    
    cpdef dict parse_attributes(self, str attr_string, str attr_format):
        """Converts a GTF or GFF3 formatted string to """
        cdef dict attr_dict = {}
        if attr_format == 'GTF':
            try:
                attr_dict = {k:v for k,v in [attr.rstrip('";').split(' "') for attr in attr_string.split('; ')]}
            except: # Failure case: attribute values aren't all surrounded by quotes
                attr_dict = {k:v.strip('"') for k,v in [attr.rstrip(';').split(' ') for attr in attr_string.split('; ')]}
        elif attr_format == 'GFF':
            attr_dict = {k:v for k,v in [attr.rstrip(';').split('=') for attr in attr_string.split(';')]}
        
        return attr_dict
    
    def __eq__(self, other): return self.span == other.span
    def __ne__(self, other): return self.span != other.span
    def __gt__(self, other): return self.span >  other.span
    def __ge__(self, other): return self.span >= other.span
    def __lt__(self, other): return self.span <  other.span
    def __le__(self, other): return self.span <= other.span


cdef class AnnotationDataset(RNAseqDataset):
    """Child of an RNAseqDataset. Parses a set of input annotations (GFF3/GTF/BED12/ELR).
    Has methods that allow these annotations to be merged together to
    form a single consensus or to be merged directly into a 'reference' annotation."""
    cdef public dict annotations, gtf_config, gff_config
    cdef public object generator
    cdef public int number_of_assemblies, counter, min_reps, confidence
    cdef public float cap_bonus
    cdef public bint verbose
    def __init__(self, annotation_files, reference=None, genome_fasta=None, config=config_defaults, gtf_config=gtf_defaults, gff_config=gff_defaults, confidence=1):
        RNAseqDataset.__init__(self, None, None, None, genome_fasta, config)
        self.min_reps = config['min_reps']
        self.cap_bonus = config['cap_bonus']
        self.verbose = config.get('verbose',False)
        self.number_of_assemblies = len(annotation_files)
        self.confidence = confidence
        self.gtf_config = gtf_config
        self.gff_config = gff_config
        cdef str f, name, k
        cdef RNAseqMapping v
        self.annotations = {}
        for f in annotation_files:
            name = f.lower().replace('.gff','').replace('.gff3','').replace('.gtf','').split('/')[-1]
            self.annotations[name] = self.import_annotation(f, name)
        
        if reference is not None:
            self.annotations['reference'] = self.import_annotation(reference, 'reference')
            for k in self.annotations['reference'].keys():
                for v in self.annotations['reference'][k]:
                    v.is_reference = True
                    v.attributes['TPM'] = 1
        
        self.counter = 0
        self.generator = self.generate_loci()

    cpdef dict import_annotation(self, str filename, str name):
        """Given a file path to a valid GTF/GFF3/BED/ELR file,
        Converts the entire file to a dict of RNAseqMapping objects.
        Each key:value pair is a chromosome:position-sorted list of objects."""
        cdef:
            RNAseqMapping item
            AnnotationObject line_object
            str file_extension, format, chrom
            dict object_dict, config_dict
            AnnotationObject current_object, current_parent
            list children
            int counter
            float total_coverage, total_s, total_e
        
        if self.verbose:
            print('Importing {}...'.format(filename))
        
        total_coverage = 0
        total_s = 0
        total_e = 0
        object_dict = {k:[] for k in self.chrom_array}
        config_dict = {}
        file_extension = filename.split('.')[-1].upper()
        if file_extension not in ['GTF','GFF3','GFF','BED','ELR']:
            return object_dict
        elif file_extension in ['GFF3','GFF']:
            format = 'GFF'
            config_dict = self.gff_config
        elif file_extension == 'GTF':
            format = 'GTF'
            config_dict = self.gtf_config
        elif file_extension == 'BED':
            format = 'BED'
        elif file_extension == 'ELR':
            format = 'ELR'
        
        file = open(filename,'r')
        current_parent = AnnotationObject('', format, config_dict) # An empty annotation object
        children = []
        if format in ['GFF','GTF']: # Parse a GFF/GTF file
            for line in file:
                if line[0] == '#':continue
                current_object = AnnotationObject(line, format, config_dict)
                if current_object.keep:
                    if current_object.parent: # Add the old object to object_dict and start a new one
                        if current_parent.parent: # Ignore the first empty annotation object
                            item = self.anno_to_mapping_object(current_parent, children, int(name=='reference'))
                            item.attributes['source'] = name
                            total_coverage += item.weight
                            total_s += float(item.attributes.get('S.reads', 0))
                            total_s += float(item.attributes.get('S.capped', 0))
                            total_e += float(item.attributes.get('E.reads', 0))
                            chrom = self.chrom_array[item.chrom]
                            if chrom not in object_dict.keys(): object_dict[chrom] = []
                            object_dict[chrom].append(item)
                        
                        current_parent = current_object
                        children = []
                    elif current_object.transcript_id == current_parent.transcript_id:
                        children.append(current_object)
            
            item = self.anno_to_mapping_object(current_parent, children, int(name=='reference'))
            item.attributes['source'] = name
            total_coverage += item.weight
            total_s += float(item.attributes.get('S.reads', 0))
            total_s += float(item.attributes.get('S.capped', 0))
            total_e += float(item.attributes.get('E.reads', 0))
            chrom = self.chrom_array[item.chrom]
            if chrom not in object_dict.keys(): object_dict[chrom] = []
            object_dict[chrom].append(item)
        elif format == 'BED':
            for line in file:
                if line[0] == '#':continue
                item = self.parse_bed_line(line)
                item.attributes['source'] = name
                total_coverage += item.weight
                chrom = self.chrom_array[item.chrom]
                if chrom not in object_dict.keys(): object_dict[chrom] = []
                object_dict[chrom].append(item)
        elif format == 'ELR':
            counter = 0
            for line in file:
                if line[0] == '#':continue
                item = self.parse_elr_line(line, name, counter)
                item.attributes['source'] = name
                total_coverage += item.weight
                chrom = self.chrom_array[item.chrom]
                if chrom not in object_dict.keys(): object_dict[chrom] = []
                object_dict[chrom].append(item)
                counter += 1
        
        file.close()
        if total_s == 0:
            total_s = total_coverage
        
        if total_e == 0:
            total_e = total_coverage
        
        for chrom in object_dict.keys():
            object_dict[chrom].sort()
            for item in object_dict[chrom]:
                item.attributes['TPM'] = item.weight/total_coverage*1000000
                if format in ['ELR','BED']:
                    item.attributes['S.reads'] = item.attributes['TPM'] if item.s_tag else 0
                    item.attributes['S.capped'] = item.attributes['TPM'] if item.capped else 0
                    item.attributes['E.reads'] = item.attributes['TPM'] if item.e_tag else 0
                
                if 'S.reads' in item.attributes:
                    item.attributes['S.ppm'] = float(item.attributes['S.reads'])/total_s*1000000
                
                if 'S.capped' in item.attributes:
                    item.attributes['C.ppm'] = float(item.attributes['S.capped'])/total_s*1000000
                
                if 'E.reads' in item.attributes:
                    item.attributes['E.ppm'] = float(item.attributes['E.reads'])/total_e*1000000
        
        return object_dict
    
    cpdef str get_transcript_fasta(self, RNAseqMapping transcript):
        """Given a transcript model, return it's mature cDNA sequence."""
        cdef:
            str chrom, fasta
        
        chrom = self.chrom_array[transcript.chrom]
        fasta = ''.join([self.genome[chrom][l:r] for l,r in transcript.ranges])
        if transcript.strand == -1:
            fasta = fu.rc(fasta)
        
        return fasta

    cdef RNAseqMapping parse_bed_line(self, str line):
        cdef RNAseqMapping new_read
        cdef list bed_elements
        input_data = parse_BED_line(line, self.chrom_dict, self.source_dict)
        if type(input_data.chrom) is str: # Chromosome wasn't in chrom_dict
            chrom_string = input_data.chrom
            input_data = input_data._replace(chrom=self.chrom_index)
            self.add_chrom(chrom_string)
        
        if type(input_data.source) is str: # Source wasn't in source_dict
            source_string = input_data.source
            input_data = input_data._replace(source=self.source_index)
            self.add_source(source_string)
        
        new_read = RNAseqMapping(input_data)
        bed_elements = line.rstrip().split('\t')
        new_read.attributes['gene_id'] = '.'.join(bed_elements[3].split('.')[:-1])
        new_read.attributes['transcript_id'] = bed_elements[3]
        new_read.attributes['S.reads'] = new_read.weight if new_read.s_tag else 0
        new_read.attributes['S.capped'] = new_read.weight if new_read.capped else 0
        new_read.attributes['E.reads'] = new_read.weight if new_read.e_tag else 0
        new_read.attributes['cov'] = new_read.weight
        return new_read
    
    cdef RNAseqMapping parse_elr_line(self, str line, str name, str counter):
        cdef RNAseqMapping new_read
        new_read = elr_to_readobject(line)
        new_read.attributes['gene_id'] = name
        new_read.attributes['transcript_id'] = '{}.{}'.format(name, counter)
        return new_read

    cpdef RNAseqMapping anno_to_mapping_object(self, AnnotationObject parent, list children, int source):
        """Given a top-level 'transcript' GTF/GFF feature and a list of 'exon' children,
        return a matching RNAseqMapping object."""
        cdef AnnotationObject child
        cdef RNAseqMapping mapping_object
        cdef list ranges, splice
        cdef int chrom, strand
        
        strand = parent.strand
        if parent.chrom not in self.chrom_dict.keys():
            self.add_chrom(parent.chrom)
        
        chrom = self.chrom_dict[parent.chrom]
        ranges = []
        children.sort()
        splice = [True]*(len(children)-1)
        for child in children:
            ranges += [child.span]

        input_data = ELdata(chrom, source, strand, ranges, splice, True, True, False, 1)
        mapping_object = RNAseqMapping(input_data, parent.attributes)
        if 'cov' in mapping_object.attributes.keys():
            mapping_object.weight = float(mapping_object.attributes['cov'])
        
        return mapping_object
    
    def generate_loci(self):
        """Yields a contiguous chunk of all annotation objects as a list
        by interleaving all sorted annotations together."""
        cdef str chrom
        cdef list current_chunk
        cdef RNAseqMapping item
        cdef int number_of_items, i, left
        number_of_files = len(self.annotations)
        for chrom in self.chrom_array:
            chrom_annotation_list = []
            for k in self.annotations.keys():
                chrom_annotation_list += self.annotations[k].get(chrom, [])
            
            chrom_annotation_list.sort()
            number_of_items = len(chrom_annotation_list)
            left = 0
            rightmost = -1
            for i in range(number_of_items):
                item = chrom_annotation_list[i]
                if rightmost > 0:
                    if item.left() > rightmost: # A gap was jumped
                        yield chrom_annotation_list[left:i]
                        left = i

                if item.span[1] > rightmost:
                    rightmost = item.span[1]
                
            yield chrom_annotation_list[left:]
        
#########################################
# Utilities for BED/ELR file processing #
#########################################

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

cpdef list get_sources(list list_of_reads):
    cdef:
        RNAseqMapping read
        set sources
        dict source_lookup
    
    sources = set([read.source for read in list_of_reads])
    return sorted(list(sources))

cpdef dict get_source_dict(list list_of_sources):
    cdef:
        set sources
        dict source_lookup
    
    sources = set(list_of_sources)
    source_lookup = dict(zip(sources, range(len(sources))))
    return source_lookup

cpdef build_depth_matrix(int leftmost, int rightmost, tuple reads, float cap_bonus=1, bint use_attributes=False):
    """Stores a numpy array of feature-specific coverage depth for the read list.
    Populates an 11-row matrix:
    S+  E+  D+  A+  S-  E-  D-  A-  cov+  cov-  cov?
    Additionally, records splice junction D-A pairs in dicts J_plus and J_minus"""
    cdef:
        Py_ssize_t Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn, covrow
        int array_length, l, r, pos
        float weight, s_weight, e_weight, c_weight
        (int, int) span, block
        str junction_hash
        dict J_plus, J_minus
        RNAseqMapping read
        np.ndarray depth_matrix
    
    Sp, Ep, Sm, Em, Dp, Ap, Dm, Am, covp, covm, covn = range(11)
    array_length = rightmost - leftmost
    depth_matrix = np.zeros(shape=(11, array_length), dtype=np.float32)
    J_plus, J_minus = {}, {}
    for read in reads:
        if use_attributes: # Check a read's attributes for different values of each type
            weight = read.weight
            s_weight = float(read.attributes.get('S.ppm', weight))
            e_weight = float(read.attributes.get('E.ppm', weight))
            c_weight = float(read.attributes.get('C.ppm', weight))
        else:
            weight = s_weight = e_weight = c_weight = read.weight
        
        if read.strand == 1:
            covrow = covp
            for span in read.junctions():
                l = span[0] - leftmost
                r = span[1] - leftmost
                depth_matrix[Dp, l] += weight
                depth_matrix[Ap, r] += weight
                block = (l,r)
                junction_hash = span_to_string(block)
                J_plus[junction_hash] = J_plus.get(junction_hash, 0) + weight
            
            if read.s_tag:
                pos = read.span[0] - leftmost
                if read.capped:
                    depth_matrix[Sp, pos] += c_weight * cap_bonus
                else:
                    depth_matrix[Sp, pos] += s_weight
            
            if read.e_tag:
                pos = read.span[1] - leftmost - 1
                depth_matrix[Ep, pos] += e_weight
        elif read.strand == -1:
            covrow = covm
            for span in read.junctions():
                l = span[0] - leftmost
                r = span[1] - leftmost
                depth_matrix[Am, l] += weight
                depth_matrix[Dm, r] += weight
                block = (l,r)
                junction_hash = span_to_string(block)
                J_minus[junction_hash] = J_minus.get(junction_hash, 0) + weight
            
            if read.e_tag:
                pos = read.span[0] - leftmost
                depth_matrix[Em, pos] += e_weight
            
            if read.s_tag:
                pos = read.span[1] - leftmost - 1
                if read.capped:
                    depth_matrix[Sm, pos] += c_weight * cap_bonus
                else:
                    depth_matrix[Sm, pos] += s_weight
        else: # The read has no features other than non-stranded coverage
            covrow = covn
        
        for span in read.ranges:
            l = span[0] - leftmost
            r = span[1] - leftmost
            depth_matrix[covrow, l:r] += weight
        
    return depth_matrix, J_plus, J_minus

cdef str span_to_string((int, int) span):
    """Converts a tuple of two ints to a string connected by ':'"""
    return '{}:{}'.format(span[0], span[1])


cdef (int, int) string_to_span(str string):
    """Converts a string from span_to_string() back into a span"""
    cdef list splitstring = string.split(':')
    return (int(splitstring[0]), int(splitstring[1]))


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
    cdef RNAseqMapping output_object = RNAseqMapping(input_data)
    return output_object



#####################################
# Utilities for BAM file processing #
#####################################
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

cpdef str get_flank(dict genome, str chrom, int pos, int strand, str label_type, int label_len):
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


cdef bint is_homopolymer(str string, float threshold=0.8):
    """Returns whether a single character composes > threshold
    of a string."""
    cdef str n
    cdef int count_n, total_count, string_length, thresh_length
    
    string = string.upper()
    string_length = len(string)
    thresh_length = int(round(string_length * threshold))
    if string_length == 0:
        return True
    
    total_count = 0
    for n in ['A','T','G','C']:
        count_n = string.count(n)
        if count_n >= thresh_length:
            return True
        else:
            total_count += count_n
            if total_count > string_length - thresh_length:
                # Enough subthreshold nucleotides were found
                return False
    
    return False


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

cdef class BAMobject:
    cdef readonly RNAseqDataset dataset
    cdef readonly list input_lines
    cdef readonly bint ignore_ends, secondary
    def __init__(self, RNAseqDataset dataset, list input_lines, bint ignore_ends=False, bint secondary=False):
        """Convert a list of pysam.AlignedSegment objects into RNAseqMappings that can be added to an RNAseqDataset.
        Quality control of end labels:
        1) An improperly mapped 5'/3' end of a read should be stripped of its tag
        2) Upstream untemplated Gs on a 5' end are evidence for a cap structure (requires genome)
        3) End labels that match the genome-templated sequence are false positive trims (requires genome)
        4) False positive oligo-dT priming can occur at genome-templated purine-rich sites (requires genome)
        5) False positive template-switching can occur at matching RNA sites of 3+ nucleotides (requires genome)
        """
        self.dataset = dataset
        self.input_lines = input_lines
        self.ignore_ends = ignore_ends
        self.secondary = secondary
    
    cpdef list generate_read(self):
        cdef:
            int s_len, s_tag_len, e_len, e_tag_len, Nmap, counter, input_len, map_number, mate, strand, number_of_blocks
            int i, gap_len, pos, junction_strand, chrom_id, start_pos, end_pos, trim_pos, errors
            float weight
            str ID, chrom, js, seq, aligned_seq, trimmed_nuc
            (bint, bint, int, int) ID_tags = (False, False, 0, 0)
            dict mappings
            list splice, gaps, ranges, introns
            bint stranded, stranded_method, fiveprime, threeprime, junction_exists
            (int, int) g
            array.array flankmatch
            RNAseqMapping current_mapping, mate_read
        
        s_tag_len = e_tag_len = 0
        stranded = stranded_method = False
        if self.dataset.stranded: # The read is strand-specific
            stranded_method = True # The method to generate the read is inherently stranded
        
        capped = self.dataset.capped
        input_len = len(self.input_lines)
        weight = float(1)/input_len
        mappings = {} # Make an empty dict to store each mapping object
        s_len = len(self.dataset.start_array)
        e_len = len(self.dataset.end_array)
        Nmap = self.get_mapping_number()
        map_number = 0
        counter = 0
        seq = ''
        if input_len == 0:
            return []

        ID = self.input_lines[0].query_name
        if not self.ignore_ends:
            ID_tags = parse_tag(ID)
        
        s_tag_len = ID_tags[2]
        e_tag_len = ID_tags[3]
        if s_tag_len > 0:
            if s_tag_len < s_len:
                s_len = s_tag_len
        
        if e_tag_len > 0:
            if e_tag_len < e_len:
                e_len = e_tag_len
        
        for i in range(input_len): 
            s_tag = self.dataset.s_tag or ID_tags[0]
            e_tag = self.dataset.e_tag or ID_tags[1]
            if s_tag or e_tag or stranded_method:
                stranded = True
            
            line = self.input_lines[i] # Each line must be a pysam.AlignedSegment
            if self.should_skip(line): # Line must pass filters
                continue
            
            mate, strand = self.determine_strand(line, stranded)
            if mate == 1:
                counter += 1
            
            try:
                map_number = line.get_tag('HI')
            except KeyError:
                map_number = counter
            
            if seq == '':
                seq = line.query_sequence
            
            pos = line.reference_start
            chrom_id = line.reference_id
            try:
                chrom = line.header.get_reference_name(chrom_id)
            except:
                chrom = line.reference_name
            
            # Parse the SAM CIGAR string to get mapped positions, splice junction sites, and softclipped positions
            try:
                errors = line.get_tag('NM')
            except KeyError:
                errors = 0
            
            try:
                mdstring = line.get_tag('MD')
            except KeyError:
                mdstring = str(len(seq))
            
            ranges, introns, head, tail = parse_SAM_CIGAR(pos, line.cigartuples, mdstring)
            number_of_blocks = len(ranges)
            if number_of_blocks == 0: # No exons of passing quality were found
                continue
            elif number_of_blocks == 1:
                alignment_strand = 0
                splice = []
            else:
                alignment_strand = self.get_alignment_strand(line)
                splice = self.get_splice_info(ranges, introns, chrom, alignment_strand) # Check which gaps between exon blocks are present in intron blocks
            
            if tail == 0:
                aligned_seq = seq[head:]
            else:
                aligned_seq = seq[head:-tail]
            
            if is_homopolymer(aligned_seq): # Aligned sequence >80% repeat of one nucleotide
                continue
            
            match_length = len(aligned_seq) - errors
            if match_length < self.dataset.minlen_loose: # Read is short enought to require stringent filtering
                if self.fails_stringent_filters(Nmap, match_length, head, tail, errors):
                    continue
            
            
            # EVALUATE SOFTCLIPPED NUCLEOTIDES
            fiveprime = mate == 1
            threeprime = (mate == 1 and not line.is_paired) or mate == 2
            s_tag, e_tag, capped = self.filter_labels_by_softclip_length(s_tag, e_tag, capped, fiveprime, threeprime, strand, head, tail)
            # Check for uuG's (5') or terminal mismatches (3')
            if self.dataset.genome:
                start_pos = 0
                end_pos = 0
                if strand == 1:
                    start_pos = ranges[0][0]
                    end_pos = ranges[-1][-1]-1
                elif strand == -1:
                    start_pos = ranges[-1][-1]-1
                    end_pos = ranges[0][0]
                
                if s_tag and fiveprime:
                    if head > 0 or tail > 0:
                        capped = self.untemplated_upstream_g(strand, head, tail, seq, chrom, ranges)
                    
                    if self.matches_masking_sequence(chrom, start_pos, strand, 'S', s_len): # Too similar to the 5' masking sequence
                        s_tag = capped = False
                
                if e_tag and threeprime:
                    if head > 0 or tail > 0:
                        self.restore_terminal_mismatches(strand, head, tail, ranges)
                    
                    if self.matches_masking_sequence(chrom, end_pos, strand, 'E', e_len): # Too similar to the 3' masking sequence
                        e_tag = False
            
            # Reconcile strand information given by start, end, and splice
            junction_exists = sum(splice) > 0
            if alignment_strand != 0 and junction_exists: # At least one strand-informative splice junction exists
                if strand != alignment_strand: # Splice disagrees with end tags; remove tags
                    strand = alignment_strand
                    s_tag = e_tag = capped = False
            
            if not stranded_method and not s_tag and not e_tag and not junction_exists:
                strand = 0 # No strand information can be found
            
            # Generate a ReadObject with the parsed attributes above
            read_data = ELdata(chrom_id, 0, strand, ranges, splice, s_tag, e_tag, capped, round(weight,2))
            current_mapping = RNAseqMapping(read_data)
            current_mapping.e_len = e_tag_len
            current_mapping.s_len = s_tag_len
            if map_number not in mappings:
                mappings[map_number] = current_mapping
            else: # merge two mate-pair ReadObjects together
                mate_read = mappings[map_number]
                mate_read.merge(current_mapping)
        
        return list(mappings.values())

    cdef int get_mapping_number(self):  
        """Given a list of pysam objects, determine
        how many locations in the genome the read (pair) mapped."""
        cdef int Nmap = 0
        cdef int num_lines = len(self.input_lines)
        if num_lines > 0: # Process the line(s)
            if num_lines == 1:
                Nmap = 1
            else: # More than one line, could be multimapper and/or paired
                line = self.input_lines[0]
                try:
                    Nmap = line.get_tag('NH')
                except KeyError:
                    if line.is_paired:
                        Nmap = int(num_lines*0.5)
                    else:
                        Nmap = num_lines
        
        return Nmap
    
    cdef (int, int) determine_strand(self, line, bint stranded):
        """Determine the RNA strand of a pysam object"""
        cdef int mate, strand
        mate = 1
        strand = 0
        if line.is_paired:
            mate = 1 if line.is_read1 else 2
        
        if mate == 1: # If stranded, sense w.r.t. RNA
            if stranded:
                strand = -1 if line.is_reverse else 1
        else:
            if stranded:
                strand = 1 if line.is_reverse else -1
        
        return mate, strand
    
    cdef bint should_skip(self, line):
        """The read should not be processed."""
        if line.is_unmapped or line.is_supplementary: # Skip unmapped reads and poorly mapped reads
            return True
        elif line.is_secondary and not self.secondary: # Unless secondary alignments are allowed, skip these too
            return True
        elif line.is_paired and not line.is_proper_pair: # Ignore discordant reads
            return True
        
        return False

    cdef bint fails_stringent_filters(self, int Nmap, int match_length, int head, int tail, int errors):
        """Reads below the 'minlen_loose' length should be treated
        more stringently: no allowed softclipping, multimapping, or mismatches.
        Absolutely require the length to be longer than minlen_strict."""
        if head > 0 or tail > 0 or errors > 0 or Nmap > 1 or match_length < self.dataset.minlen_strict:
            return True
        
        return False

    cdef int get_alignment_strand(self, line):
        """Returns 1(+), -1(-), or 0(.) if one of the BAM
        splice tags (XS, ts) contains strand information."""
        cdef str js
        cdef alignment_strand = 0
        try:
            js = line.get_tag('XS')
        except KeyError:
            try:
                js = line.get_tag('ts')
            except KeyError:
                js = '.'
        
        if js == '+':
            alignment_strand = 1
        elif js == '-':
            alignment_strand = -1
        
        return alignment_strand

    cdef list get_splice_info(self, list ranges, list introns, str chrom, int alignment_strand):
        """Returns a list of booleans denoting whether each gap
        between ranges is a splice junction or not."""
        cdef list splice
        cdef Py_ssize_t i, range_len, gap_len
        cdef junction_strand
        cdef str js = '.'
        range_len = len(ranges)
        gap_len = range_len - 1
        gaps = [(ranges[i][1], ranges[i+1][0]) for i in range(range_len-1)] # List of all gaps between ranges
        splice = [False]*gap_len # List of booleans indicating whether each gap is a splice junction
        for i in range(gap_len):
            g = gaps[i]
            if g in introns:
                splice[i] = True
                if self.dataset.genome:
                    junction_strand = get_junction_strand(self.dataset.genome, chrom, g[0], g[1])
                    if junction_strand != alignment_strand: # Did not find a valid splice junction
                        splice[i] = False
        
        return splice
    
    cdef (bint, bint, bint) filter_labels_by_softclip_length(self, bint s_tag, bint e_tag, bint capped, bint fiveprime, bint threeprime, int strand, int head, int tail):
        """Determines whether the s_tag, e_tag and capped parameters
        should be removed an alignment that has softclipping on its edges."""  
        if strand == 0:
            return False, False, False

        if not fiveprime:
            s_tag = capped = False
        else:
            if strand == 1 and (head == -1 or head > 4):
                s_tag = capped = False
            elif strand == -1 and (tail == -1 or tail > 4):
                s_tag = capped = False
        
        if not threeprime:
            e_tag = False
        else:
            if strand == 1 and (tail == -1 or tail > 4):
                e_tag = False
            elif strand == -1 and (head == -1 or head > 4):
                e_tag = False
        
        return s_tag, e_tag, capped
    
    cdef bint untemplated_upstream_g(self, int strand, int head, int tail, str seq, str chrom, list ranges):
        """Checks (1) if a softclipped string at a read's 5' end
        is an oligomer of G and (2) if that oligomer does not match the genome.
        Returns True if evidence supports a cap."""
        if strand == 1:
            if head <= 0 or head > 4:
                return False
            
            if seq[:head] == 'G'*head: # Softclipped nucleotides are G
                if get_flank(self.dataset.genome, chrom, ranges[0][0], 1, 'S', 1) != 'G': # The flanking nucleotide is NOT G
                    return True # One or more upstream untemplated Gs were detected
                else:
                    return False
        elif strand == -1:
            if tail <= 0 or tail > 4:
                return False
            
            if seq[-tail:] == 'C'*tail: # Sofclipped nucleotides are (antisense) G
                if get_flank(self.dataset.genome, chrom, ranges[-1][-1]-1, -1, 'S', 1) != 'G': # The flanking nucleotide is NOT G
                    return True
                else:
                    return False

    cdef void restore_terminal_mismatches(self, int strand, int head, int tail, list ranges):
        """Updates the mapping ranges of a read with a softclipped
        sequenced added back to one end."""
        if strand == 1: 
            if tail > 0: # Right clip exists
                ranges[-1] = (ranges[-1][0],ranges[-1][1]+tail)
        elif strand == -1:
            if head > 0: # Left clip exists
                ranges[0] = (ranges[0][0]-head,ranges[0][1])

    cdef bint matches_masking_sequence(self, str chrom, int position, int strand, str readtype, int length):
        """Evaluates whether a clipped tag matches too closely
        with a genome-templated region could have caused 
        false positive end signal"""
        ## 5'
        flank = get_flank(self.dataset.genome, chrom, position, strand, readtype, length) # Get upstream flanking sequence to start
        if len(flank) > 0:
            if readtype == 'S':
                flankmatch = self.dataset.start_array[-length:]
            elif readtype == 'E':
                flankmatch = self.dataset.end_array[:length]
                
            if fu.IUPACham(fu.nuc_to_int(flank), flankmatch, self.dataset.mismatch_rate*length) <= self.dataset.mismatch_rate*length:
                return True
        
        return False


def read_generator(fileconn, RNAseqDataset dataset, str file_type, int max_gap, float minimum_proportion):
    """Yields a contiguous chunk of reads from the input file
    separated on either side by a gaps > max_gap"""
    cdef RNAseqMapping last_read
    cdef int l, r, old_chrom, old_l, old_r, rightmost, k
    cdef (int, int) exon
    cdef float read_bases, current_bases, mean_cov, current_cov, total_cov
    cdef set covered_positions
    cdef list passed_positions
    if file_type == 'elr':
        add_read = dataset.add_read_from_ELR
    elif file_type == 'bed':
        add_read = dataset.add_read_from_BED
    elif file_type == 'bam' or file_type == 'sam':
        add_read = dataset.add_read_from_BAM
    else:
        return
    
    end_positions = Counter() # Keep track of where reads end to maintain a tally of coverage depth
    old_chrom, old_l, old_r, rightmost = -1, -1, -1, -1
    current_cov = 0
    current_bases = 0
    mean_cov = 0
    covered_positions = set()
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
        read = dataset.read_list[-1]
        read_bases = read.get_length()*read.weight
        l, r = read.span
        current_cov += read.weight
        if old_chrom == -1: # Uninitialized; add the read and make no other decisions
            pass
        elif read.chrom != old_chrom or l >= rightmost + max_gap: # The last locus is definitely finished; dump the read list
            yield dataset.read_list[:-1]
            dataset.read_list = [read]
            covered_positions = set()
            for exon in read.ranges:
                covered_positions.update(range(exon[0],exon[1]))
            
            current_bases = 0
            current_cov = read.weight
            mean_cov = current_cov
            end_positions = Counter()
            rightmost = r
        elif l > old_l: # Read advanced, but not by enough to automatically cut
            passed_positions = [k for k in end_positions.keys() if k <= l]
            for k in passed_positions:
                current_cov -= end_positions.pop(k)
            
            if current_cov < minimum_proportion * mean_cov: # Current cov is sufficiently low to cause a break
                yield dataset.read_list[:-1]
                dataset.read_list = [read]
                covered_positions = set()
                for exon in read.ranges:
                    covered_positions.update(range(exon[0],exon[1]))
                
                current_bases = 0
                current_cov = read.weight
                end_positions = Counter()
                rightmost = r
        
        end_positions[r] += read.weight # Add the read's weight to the position where the read ends
        current_bases += read_bases
        mean_cov = current_bases / len(covered_positions)
        if r > rightmost: rightmost = r
        old_chrom, old_l, old_r = read.chrom, l, r
    
    # Dump the remaining reads
    yield dataset.read_list
    fileconn.close()


cpdef (int, int) get_max_deltas(np.ndarray[float, ndim=1] array, float offset):
    """Returns two positions in the array that represent
      [0] the sharpest increase and 
      [1] the sharpest decrease
    in the float values in the array.
    """
    cdef Py_ssize_t i, j, l
    cdef float vi, vj, delta, max_delta, min_delta
    cdef (int, int) top_positions = (0, 0)
    cdef float [:] ARRAY = array
    if offset < 1:
        offset = 1
    
    l = array.shape[0]
    max_delta = min_delta = 0
    for i in range(1, l):
        j = i - 1
        vi = ARRAY[i]
        vj = ARRAY[j]
        delta = (vi - vj)/(vi + offset)
        if delta > max_delta: # A gap begins or extends
            max_delta = delta
            top_positions[0] = i
        elif delta < min_delta: # Not in a gap
            min_delta = delta
            top_positions[1] = i-1
    
    return top_positions


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

