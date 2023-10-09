#cython: language_level=3
cimport cython
import array
import sys
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
ELdata = namedtuple('ELdata', 'chrom source strand ranges splice s_tag e_tag capped weight condensed')
cdef class RNAseqMapping():
    cdef public int chrom, source, strand, s_len, e_len
    cdef public list ranges, splice
    cdef public bint s_tag, e_tag, capped, complete, is_reference, condensed
    cdef public dict attributes
    cdef public (int, int) span
    cdef public float weight, coverage
    cdef public (float, float, float) tripleweight
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
        self.condensed = input_data.condensed
        if self.strand == 0: # Terminal tag information is meaningless for nonstranded reads
            self.s_tag = self.e_tag = self.capped = False
        else:
            self.s_tag, self.e_tag, self.capped = input_data.s_tag, input_data.e_tag, input_data.capped
        
        self.span = (self.left(), self.right())
        self.complete = False
        if self.s_tag and self.e_tag and False not in self.splice:
            self.complete = True
        
        if attributes is not None:
            self.attributes = attributes
        else:
            self.attributes = {}
        
        if '|' in str(input_data.weight):
            tripleweight = tuple(float(s) for s in str(input_data.weight).split('|'))
            self.weight = tripleweight[0]
            self.attributes['E.reads'] = tripleweight[2]
            if self.capped:
                self.attributes['S.capped'] = tripleweight[1]
            else:
                self.attributes['S.reads'] = tripleweight[1]
        else:
            self.weight = float(input_data.weight)
    
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
        return (max([self.span[0],other.span[0]]), min([self.span[1],other.span[1]]))
    
    cpdef int shared_bases(self, RNAseqMapping other):
        """Returns the number of exonic nucleotides that overlap between
        self and other."""
        cdef int index_self, index_other, exlen_self, exlen_other, l1, r1, l2, r2, shared
        if not self.overlaps(other):
            return 0
        
        shared = 0
        exlen_self = len(self.ranges)
        exlen_other = len(other.ranges)
        index_self, index_other = 0, 0
        while index_self < exlen_self and index_other < exlen_other:
            l1, r1 = self.ranges[index_self]
            l2, r2 = other.ranges[index_other]
            if l2 > r1: # other is fully right of self
                index_self += 1
                continue
            
            if l1 > r2: # self is fully right of other
                index_other += 1
                continue
            
            shared += max(0, min(r1,r2)-max(l1,l2))
            if r1 > r2: # self is right of other, advance other
                index_other += 1
            elif r2 > r1: # other is right of self, advance self
                index_self += 1
            else: # right edges are the same, advance both
                index_other += 1
                index_self += 1
        
        return shared
    
    cpdef bint is_malformed(self):
        """Returns if any exons are length zero."""
        cdef int l,r
        for l,r in self.ranges:
            if r <= l:
                return True
        
        return False
    
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

    cpdef bint is_compatible(self, RNAseqMapping other, bint ignore_ends=False, bint ignore_source=False, bint ignore_overhang=False, bint allow_intron_retention=False):
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
        if ignore_overhang:
            j1 = [j for j in self.junctions() if j[0] > overlap[0] and j[1] < overlap[1]]
            j2 = [j for j in other.junctions() if j[0] > overlap[0] and j[1] < overlap[1]]
        else:
            j1 = [j for j in self.junctions() if j[1] > overlap[0] and j[0] < overlap[1]]
            j2 = [j for j in other.junctions() if j[1] > overlap[0] and j[0] < overlap[1]]
        
        if j1 == j2:
            return True
        
        if allow_intron_retention and all([j in j2 for j in j1]):
            # Every intron contained in j1 is also present in j2
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

    cpdef str get_node_labels(self, bint record_artifacts=False, bint condense=False):
        """Returns a string with one label for each edge of each range in self.ranges."""
        cdef str startchar, endchar, gapchar
        gapchar = ['AD','..','DA'][1+self.strand]
        startchar = ['.',['S','C'][self.capped]][int(self.s_tag or (record_artifacts and self.s_len > 0))]
        endchar = ['.','E'][int(self.e_tag or (record_artifacts and self.e_len > 0))]
        if record_artifacts:
            if startchar != '.' and not self.s_tag:
                startchar = '>'
            
            if endchar != '.'  and not self.e_tag:
                endchar = '<'
        
        if self.strand == -1:
            startchar, endchar = endchar, startchar

        if condense:
            startchar = startchar.lower()
            endchar = endchar.lower()
        
        return ''.join([startchar]+[gapchar if i else '..' for i in self.splice]+[endchar])
    
    cpdef write_as_elr(self, bint as_string=True, bint record_artifacts=False, bint condense=False, bint endweights=False):
        """Returns a string that represents the ReadObject
        in the end-labeled read (ELR) format"""
        cdef str elr_strand, labels, weightstring
        cdef list block_ends, elr_line
        elr_strand = '.'
        if self.strand == 1:
            elr_strand = '+'
        elif self.strand == -1:
            elr_strand = '-'
        
        block_ends = flatten(self.ranges)
        lengths = [block_ends[i]-block_ends[i-1] for i in range(1,len(block_ends))]
        labels = self.get_node_labels(record_artifacts, condense)
        EL_CIGAR = ''.join([str(a)+str(b) for a,b in zip(labels,lengths+[''])])
        read_len = self.right() - self.left()
        if endweights and (self.s_tag or self.e_tag):
            weightstring = '{}|{}|{}'.format(
                round(self.weight,2),
                round(float(self.attributes.get('S.reads', 0))+float(self.attributes.get('S.capped', 0)),2),
                round(float(self.attributes.get('E.reads', 0)),2)
            )
        else:
            weightstring = str(round(self.weight,2))
        
        elr_line = [self.chrom, self.left(), read_len, elr_strand, EL_CIGAR, self.source, weightstring]
        if as_string:
            return '\t'.join([str(i) for i in elr_line])
        else:
            return elr_line
     
    cpdef write_as_bed(self, chrom_array, source_array, as_string=True, score_column='weight', record_artifacts=False, name_attr=None, color=None, condense=False, longStart=None, longEnd=None):
        """Returns a string that represents the ReadObject
        in a 15-column BED format"""
        cdef int chromStart, chromEnd
        labels = self.get_node_labels(record_artifacts, condense)
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
            elif ends in ['CE','EC']:
                rgb = bed_colors['CE']
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
        
        chromStart = self.ranges[0][0]
        chromEnd = self.right()
        if longStart is None:
            longStart = chromStart
        else:
            self.ranges[0] = (min(self.ranges[0][0],longStart), self.ranges[0][1])
        
        if longEnd is None:
            longEnd = chromEnd
        else:
            self.ranges[-1] = (self.ranges[-1][0], max(self.ranges[-1][1],longEnd))
        
        longStart, blockStarts, blockSizes = explode_block_ranges(self.ranges)
        bed_line = [
            chrom_array[self.chrom], longStart, longEnd,
            name, score, bed_strand, chromStart, chromEnd, rgb,
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
    'reverse':False,
    'start_seq':'ACGGG',
    'end_seq':'RRRRRRRRRRRRRRRRRRRR',
    'minlen_strict':20,
    'minlen_loose':25,
    'mismatch_rate':0.2,
    'error_rate':0.1,
    'min_reps':2,
    'cap_bonus':5,
    'sj_shift':0,
    'max_headclip':4,
    'max_intron':100000,
    'confidence_threshold':0.5,
    'gene_delim':'.',
    'remove_noncanonical':False,
    'remove_gapped_termini':False,
    'secondary':False,
    'ignore_ends':False,
    'labels_are_trimmed':True,
    'quality_filter':True,
    'reference':None,
    'sj':None,
}
gtf_defaults = {
    'parent_types':set(['transcript']),
    'parent_key_transcript':['transcript_id', 'Name'],
    'parent_key_gene':'gene_id',
    'child_types':set(['exon']),
    'child_key_transcript':['transcript_id', 'Parent'],
    'child_key_gene':'gene_id'
}
gff_defaults = {
    'parent_types':set([    
        'mRNA','transcript',
        'snoRNA','tRNA','snRNA', 'miRNA','rRNA','ncRNA','mRNA_TE_gene','pseudogenic_transcript',
        'antisense_lncRNA','antisense_RNA','lnc_RNA', 'primary_transcript',
        'guide_RNA', 'scRNA', 'RNase_MRP_RNA', 'Y_RNA', 'RNase_P_RNA', 'telomerase_RNA']),
    'parent_key_transcript':['transcript_id'],
    'parent_key_gene':'gene',
    'child_types':set(['exon','pseudogenic_exon']),
    'child_key_transcript':['transcript_id', 'Parent'],
    'child_key_gene':'gene'
}
cdef class RNAseqDataset():
    cdef public list read_list, chrom_array, source_array, chrom_lengths
    cdef public int chrom_index, source_index, sj_shift, max_headclip, max_intron
    cdef public dict chrom_dict, source_dict, reference_dict, gtf_config, gff_config
    cdef readonly dict config, genome, label_tally
    cdef readonly bint s_tag, e_tag, capped, stranded, reverse, ignore_ends, remove_noncanonical, labels_are_trimmed, quality_filter, verbose, secondary, has_genome, remove_gapped_termini
    cdef readonly str start_seq, end_seq, start_seq_rc, end_seq_rc, gene_delim
    cdef readonly int minlen, minlen_strict, minlen_loose
    cdef readonly float mismatch_rate, error_rate
    cdef readonly array.array start_array, end_array
    cdef public set sj_set

    def __init__(self, chrom_array=None, source_array=None, chrom_lengths=None, genome_fasta=None, config=config_defaults, gtf_config=gtf_defaults, gff_config=gff_defaults):
        """Container for RNAseqMapping objects. Stores a reference dictionary for all
        chromosome names and sample names. Contains methods for parsing
        a variety of files into a collection of read objects."""
        self.config = config_defaults
        self.config.update(config)
        self.label_tally = {'S':Counter(), 's':Counter(), 'E':Counter(), 'e':Counter()}
        self.read_list = []
        self.gtf_config = gtf_config
        self.gff_config = gff_config
        self.s_tag = self.config['s_tag']
        self.e_tag = self.config['e_tag']
        self.capped = self.config['capped']
        self.stranded = self.config['stranded']
        self.reverse = self.config['reverse']
        self.start_seq = str(self.config['start_seq'])
        self.end_seq = str(self.config['end_seq'])
        self.minlen_strict = self.config['minlen_strict']
        self.minlen_loose = self.config['minlen_loose']
        self.minlen = self.minlen_strict
        self.max_intron = self.config['max_intron']
        self.mismatch_rate = self.config['mismatch_rate']
        self.error_rate = self.config['error_rate']
        self.sj_shift = self.config['sj_shift']
        self.max_headclip = self.config['max_headclip']
        self.remove_noncanonical = self.config['remove_noncanonical']
        self.remove_gapped_termini = self.config['remove_gapped_termini']
        self.labels_are_trimmed = self.config['labels_are_trimmed']
        self.quality_filter = self.config['quality_filter']
        self.ignore_ends = self.config['ignore_ends']
        self.secondary = self.config['secondary']
        self.gene_delim = self.config['gene_delim']
        self.verbose = self.config.get('verbose', False)
        if self.start_seq in ['None', 'False', '']:
            self.start_seq_rc = ''
            self.start_array = fu.nuc_to_int('x')
        else:    
            self.start_seq_rc = fu.rc(self.start_seq)
            self.start_array = fu.nuc_to_int(self.start_seq)
        
        if self.end_seq in ['None', 'False', '']:
            self.end_seq_rc = ''
            self.end_array = fu.nuc_to_int('x')
        else:    
            self.end_seq_rc = fu.rc(self.end_seq)
            self.end_array = fu.nuc_to_int(self.end_seq)
        
        self.chrom_lengths = chrom_lengths
        self.reference_dict = {}
        self.chrom_dict = {}
        self.chrom_index = 0
        self.chrom_array = []
        if chrom_array is not None:
            for c in chrom_array:
                self.add_chrom(c)
        
        if genome_fasta is not None:
            self.genome, index = fu.import_genome(genome_fasta)
            index_lines = [l.split('\t') for l in index.rstrip().split('\n')]
            self.chrom_array = [l[0] for l in index_lines]
            self.chrom_lengths = [int(l[1]) for l in index_lines]
            self.chrom_index = len(self.chrom_array)
            self.chrom_dict = dict(zip(self.chrom_array, range(self.chrom_index)))
            self.has_genome = True
        else:
            self.genome = {}
            self.has_genome = False
        
        self.source_dict = {}
        self.source_index = 0
        self.source_array = []
        if source_array is not None:
            for s in source_array:
                self.add_source(s)
        
        if self.config.get('reference', None) is not None:
            self.reference_dict = self.import_annotation(self.config['reference'], 'reference')
        
        self.sj_set = self.make_sj_set(self.config.get('sj',None))
    
    cpdef make_sj_set(self, sj_file):
        '''Create a reference set of known splice junctions
        with provided annotation file(s).'''
        cdef RNAseqMapping transcript
        cdef tuple sj_tuple
        cdef int l, r
        cdef list fields
        sj_set = set()
        if sj_file is not None:
            testline = open(sj_file,'r').readline()
            fields = testline.split('\t')
            if len(fields) == 9 or sj_file.lower().endswith('.out.tab'): # STAR SJ.out.tab file
                for line in open(sj_file,'r'):
                    fields = line.rstrip().split('\t')
                    if fields[3]=='1': # Forward strand
                        sj_tuple = (fields[0], int(fields[1])-1, int(fields[2]), 1)
                    elif fields[3]=='2': # Reverse strand
                        sj_tuple = (fields[0], int(fields[1])-1, int(fields[2]), -1)
                    
                    sj_set.add(sj_tuple)
            elif len(fields) == 6 and sj_file.lower().endswith('.bed'): # intron BED6 file
                for line in open(sj_file,'r'):
                    fields = line.rstrip().split('\t')
                    if fields[5]=='+': # Forward strand
                        sj_tuple = (fields[0], int(fields[1]), int(fields[2]), 1)
                    elif fields[5]=='-': # Reverse strand
                        sj_tuple = (fields[0], int(fields[1]), int(fields[2]), -1)
                    
                    sj_set.add(sj_tuple)
            else:
                print("ERROR: --splice file not recognized. Provide SJ.out.tab or BED6.")
                print("For BED12/GTF/GFF3, use --reference.")
                sys.exit(1)
        
        for chrom in self.reference_dict.keys():
            for transcript in self.reference_dict[chrom]:
                for l,r in transcript.junctions():
                    sj_tuple = (chrom, l, r, transcript.strand)
                    sj_set.add(sj_tuple)

        return sj_set

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
    
    cpdef add_read_from_BAM(self, bam_lines):
        cdef list new_read_list
        cdef RNAseqMapping read
        cdef BAMobject BAM
        if type(bam_lines) is not list:
            bam_lines = [bam_lines]
        
        BAM = BAMobject(dataset=self, input_lines=bam_lines)
        new_read_list = BAM.generate_read()
        if len(new_read_list) > 0:
            read = new_read_list[0]
            if read.s_tag:
                self.label_tally['S'][read.s_len] += 1
            elif read.s_len > 0:
                self.label_tally['s'][read.s_len] += 1
            
            if read.e_tag:
                self.label_tally['E'][read.e_len] += 1
            elif read.e_len > 0:
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
    
    cpdef dict import_annotation(self, str filename, str name):
        """Given a file path to a valid GTF/GFF3/BED/ELR file,
        Converts the entire file to a dict of RNAseqMapping objects.
        Each key:value pair is a chromosome:position-sorted list of objects."""
        cdef:
            RNAseqMapping item
            AnnotationObject current_object, current_parent, last_child
            str file_extension, format, chrom
            dict object_dict, config_dict
            list children
            int counter
            float total_coverage, total_s, total_e, coverage, s, e
        
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
        if format in ['GFF','GTF']: # Parse a GFF/GTF file
            current_parent = AnnotationObject('', format, config_dict) # An empty annotation object
            current_object = AnnotationObject('', format, config_dict) # An empty annotation object
            last_child = current_object
            children = []
            for line in file:
                if line[0] == '#':continue
                current_object = AnnotationObject(line, format, config_dict)
                if current_object.keep:
                    if current_object.parent: # Add the old object to object_dict and start a new one
                        if current_parent.parent: # Ignore the first empty annotation object
                            coverage, s, e = self.add_mapping_object(current_parent, children, name, int(name=='reference'), object_dict)
                            total_coverage += coverage
                            total_s += s
                            total_e += e
                        
                        current_parent = current_object
                        children = []
                    else:
                        if current_object.transcript_id and current_object.transcript_id == last_child.transcript_id or len(children) == 0:
                            children.append(current_object)
                        else: # New transcript from same parent
                            coverage, s, e = self.add_mapping_object(current_parent, children, name, int(name=='reference'), object_dict)
                            total_coverage += coverage
                            total_s += s
                            total_e += e
                            children = [current_object]
                        
                        last_child = current_object
            
            coverage, s, e = self.add_mapping_object(current_parent, children, name, int(name=='reference'), object_dict)
            total_coverage += coverage
            total_s += s
            total_e += e

        elif format == 'BED':
            for line in file:
                if line[0] == '#':continue
                item = self.parse_bed_line(line, name)
                item.attributes['source'] = name
                total_coverage += item.weight
                chrom = self.chrom_array[item.chrom]
                if chrom not in object_dict.keys():
                    object_dict[chrom] = []
                
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
                item.attributes['TPM'] = round(item.weight/total_coverage*1000000,2) if total_coverage>0 else 0
                if format in ['ELR','BED']:
                    item.attributes['S.reads'] = item.attributes['TPM'] if item.s_tag else 0
                    item.attributes['S.capped'] = item.attributes['TPM'] if item.capped else 0
                    item.attributes['E.reads'] = item.attributes['TPM'] if item.e_tag else 0
        
        return object_dict
    
    cpdef (float, float, float) add_mapping_object(self, AnnotationObject parent, list children, str name, int source, dict object_dict):
        """Converts a GTF/GFF collection of AnnotationObjects to a single RNAseqMapping object"""
        cdef RNAseqMapping item
        cdef AnnotationObject child
        cdef str transcript_id, chrom
        if len(children) == 0:return (0., 0., 0.)
        if parent.keep == False:return (0., 0., 0.)
        transcript_id = children[0].transcript_id
        if not transcript_id:
            return (0., 0., 0.)
        
        for child in children:
            if len(child.gene_id)>0 and child.gene_id != parent.gene_id:return (0., 0., 0.)
            if child.transcript_id != transcript_id:return (0., 0., 0.)
        
        item = self.anno_to_mapping_object(parent, children, source)
        item.attributes['source'] = name
        item.attributes['transcript_id'] = transcript_id
        chrom = self.chrom_array[item.chrom]
        if chrom not in object_dict.keys(): object_dict[chrom] = []
        object_dict[chrom].append(item)
        return (item.weight, float(item.attributes.get('S.reads', 0))+float(item.attributes.get('S.capped', 0)), float(item.attributes.get('E.reads', 0)))

    cpdef RNAseqMapping anno_to_mapping_object(self, AnnotationObject parent, list children, int source):
        """Given a top-level 'transcript' GTF/GFF feature and a list of 'exon' children,
        return a matching RNAseqMapping object."""
        cdef AnnotationObject child
        cdef RNAseqMapping mapping_object
        cdef list ranges, splice
        cdef int chrom, strand
        
        children = self.merge_children(children)
        strand = parent.strand
        if parent.chrom not in self.chrom_dict.keys():
            self.add_chrom(parent.chrom)
        
        chrom = self.chrom_dict[parent.chrom]
        ranges = []
        children.sort()
        splice = [True]*(len(children)-1)
        for child in children:
            ranges += [child.span]

        input_data = ELdata(chrom, source, strand, ranges, splice, True, True, False, 1, False)
        mapping_object = RNAseqMapping(input_data, copy.copy(parent.attributes))
        mapping_object.attributes['transcript_id'] = child.transcript_id
        if 'cov' in mapping_object.attributes.keys():
            mapping_object.weight = float(mapping_object.attributes['cov'])
        
        return mapping_object
    
    cdef RNAseqMapping parse_bed_line(self, str line, str source_string):
        cdef RNAseqMapping new_read
        cdef list bed_elements
        input_data = parse_BED_line(line, self.chrom_dict, self.source_dict, source_string, s_tag=True, e_tag=True, capped=True, gaps_are_junctions=True, keep_readname=True)
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
        new_read.attributes['gene_id'] = '.'.join(bed_elements[3].split(self.gene_delim)[:-1])
        new_read.attributes['transcript_id'] = bed_elements[3]
        new_read.attributes['S.reads'] = new_read.weight if new_read.s_tag else 0
        new_read.attributes['S.capped'] = new_read.weight if new_read.capped else 0
        new_read.attributes['E.reads'] = new_read.weight if new_read.e_tag else 0
        new_read.attributes['cov'] = new_read.weight
        return new_read
    
    cdef RNAseqMapping parse_elr_line(self, str line, str name, int counter):
        cdef RNAseqMapping new_read
        new_read = elr_to_readobject(line)
        new_read.attributes['gene_id'] = name
        new_read.attributes['transcript_id'] = '{}.{}'.format(name, counter)
        return new_read



######################################
# Utilities for GTF/GFF file parsing #
######################################
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
    'lnc_RNA': '56,114,168', 
    'macro_lncRNA': '6,74,140', 
    'miRNA': '198,95,84', 
    'misc_RNA': '249,185,54', 
    'Mt_rRNA': '0,0,0', 
    'Mt_tRNA': '0,0,0', 
    'mRNA_TE_gene': '120,120,120',
    'ncRNA': '56,114,168',
    'nonsense_mediated_decay': '180,155,100', 
    'non_stop_decay': '180,155,100', 
    'nontranslating_CDS': '180,155,100',
    'otherRNA': '56,114,168',
    'polymorphic_pseudogene': '80,80,80', 
    'pseudogenic_transcript': '80,80,80', 
    'pre_miRNA': '198,95,84', 
    'processed_pseudogene': '80,80,80', 
    'processed_transcript': '180,155,100', 
    'primary_transcript': '198,95,84',
    'protein_coding': '49,132,44', 
    'mRNA': '49,132,44', 
    'transcript': '0,0,0', 
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
            str gene_id_key, t_id
            list transcript_id_keys
        
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
                        gene_id_key = config_dict['parent_key_gene']
                        transcript_id_keys = config_dict['parent_key_transcript']
                    else:
                        gene_id_key = config_dict['child_key_gene']
                        transcript_id_keys = config_dict['child_key_transcript']
                    
                    self.gene_id = self.attributes.get(gene_id_key,'')
                    self.gene_id = self.gene_id.split(':')[-1]
                    self.transcript_id = ''
                    for t_id in transcript_id_keys:
                        self.transcript_id = self.attributes.get(t_id,'')
                        if self.transcript_id != '':
                            self.transcript_id = self.transcript_id.split(':')[-1]
                            break
                    
                    self.attributes['gene_id'] = self.gene_id
                    self.attributes['transcript_id'] = self.transcript_id
                    self.attributes['anno_type'] = self.anno_type
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
                attr_dict = {k:v.rstrip('";') for k,v in [attr.split(' "') for attr in attr_string.split('"; ')]}
            except: # Failure case: attribute values aren't all surrounded by quotes
                attr_dict = {k:v.strip('"') for k,v in [attr.rstrip(';').split(' ') for attr in attr_string.split('; ')]}
        elif attr_format == 'GFF':
            attr_dict = {k:v for k,v in [attr.rstrip(';').split('=') for attr in attr_string.rstrip(';').split(';')]}
        
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
    cdef public dict annotations
    cdef public object generator
    cdef public int number_of_assemblies, counter, min_reps, confidence
    cdef public float cap_bonus
    def __init__(self, annotation_files, reference=None, genome_fasta=None, config=config_defaults, gtf_config=gtf_defaults, gff_config=gff_defaults, confidence=1):
        RNAseqDataset.__init__(self, None, None, None, genome_fasta, config)
        self.source_array = ['reference'] + annotation_files
        self.min_reps = config['min_reps']
        self.cap_bonus = config['cap_bonus']
        self.verbose = config.get('verbose',False)
        self.gene_delim = self.config['gene_delim']
        self.number_of_assemblies = len(annotation_files)
        self.confidence = confidence
        self.gtf_config = gtf_config
        self.gff_config = gff_config
        cdef RNAseqMapping v
        self.chrom_dict = {k:i for k,i in zip(self.chrom_array,range(len(self.chrom_array)))}
        self.source_dict = {k:i for k,i in zip(self.source_array,range(len(self.source_array)))}
        self.annotations = {}
        if reference is not None:
            self.annotations['reference'] = self.import_annotation(reference, 'reference')
            for k in self.annotations['reference'].keys():
                for v in self.annotations['reference'][k]:
                    v.is_reference = True
                    v.attributes['TPM'] = 1
        
        for f in annotation_files:
            self.annotations[f] = self.import_annotation(f, f)
        
        self.counter = 0
        self.generator = self.generate_loci()

    cpdef str get_transcript_fasta(self, RNAseqMapping transcript):
        """Given a transcript model, return it's mature cDNA sequence."""
        cdef:
            str chrom, fasta
            int leftExtend, rightExtend
            list ranges
        
        ranges = copy.deepcopy(transcript.ranges)
        chrom = self.chrom_array[transcript.chrom]
        leftExtend = 0
        rightExtend = 0
        if transcript.strand == 1:
            if 'S.left' in transcript.attributes:
                leftExtend = ranges[0][0] - int(transcript.attributes['S.left'])
                ranges[0] = (ranges[0][0]-leftExtend, ranges[0][1])
            
            if 'E.right' in transcript.attributes:
                rightExtend = int(transcript.attributes['E.right']) - ranges[-1][1]
                ranges[-1] = (ranges[-1][0], ranges[-1][1]+rightExtend)
        elif transcript.strand == -1:
            if 'E.left' in transcript.attributes:
                leftExtend = ranges[0][0] - int(transcript.attributes['E.left'])
                ranges[0] = (ranges[0][0]-leftExtend, ranges[0][1])
            
            if 'S.right' in transcript.attributes:
                rightExtend = int(transcript.attributes['S.right']) - ranges[-1][1]
                ranges[-1] = (ranges[-1][0], ranges[-1][1]+rightExtend)
        
        fasta = ''.join([self.genome[chrom][l:r] for l,r in ranges])
        if leftExtend>0:
            fasta = fasta[:leftExtend].lower()+fasta[leftExtend:]

        if rightExtend>0:
            fasta = fasta[:-rightExtend]+fasta[-rightExtend:].lower()
        
        if transcript.strand == -1:
            fasta = fu.rc(fasta)
        elif transcript.strand == 0:
            fasta = fasta.lower()
        
        return fasta
    
    cpdef list merge_children(self, children):
        """Checks if there is a distance of 0 between any adjacent children;
        if so, treats them as one object."""
        cdef AnnotationObject child, last_child
        cdef list children_out
        children_out = []
        for child in sorted(children):
            if len(children_out) == 0:
                children_out.append(child)
                last_child = child
            else:
                if child.span[0] == last_child.span[1]:
                    last_child.span[1] = child.span[1]
                else:
                    children_out.append(child)
                    last_child = child
        
        return children_out
    
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
    'U':'139,137,138',
    'C':'34,82,37',
    'S':'79,153,186',
    'E':'215,60,41',
    'SE':'155,105,178',
    'CE':'25,5,23'
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
    sizes = [int(s) for s in blockSizes.rstrip(',').split(',')]
    starts = [int(s) for s in blockStarts.rstrip(',').split(',')]
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
    cdef str chrom_num, chromStart, read_len, elr_strand, EL_CIGAR, source_num, weight_string, first, last
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
    startpos = int(chromStart)
    label_indices = [i for i,character in enumerate(EL_CIGAR) if not (character.isdigit() or character=='-')]
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
    first = EL_CIGAR[0].upper()
    last = EL_CIGAR[-1].upper()
    condensed = first != EL_CIGAR[0] or last != EL_CIGAR[-1]
    s_tag = first in ['C', 'S'] or last in ['C','S']
    capped = first == 'C' or last == 'C'
    e_tag = first == 'E' or last == 'E'
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight_string, condensed)

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

cpdef build_depth_matrix(int leftmost, int rightmost, tuple reads, bint use_attributes=True, bint splice=True):
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
    
    Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
    array_length = rightmost - leftmost
    depth_matrix = np.zeros(shape=(9, array_length), dtype=np.float32)
    J_plus, J_minus = {}, {}
    for read in reads:
        if read.condensed:
            s_weight = 1 if read.s_tag else 0
            c_weight = 1 if read.capped else 0
            e_weight = 1 if read.e_tag else 0
            weight = read.weight
        elif use_attributes: # Check a read's attributes for different values of each type
            weight = read.weight
            s_weight = float(read.attributes.get('S.reads', weight))
            e_weight = float(read.attributes.get('E.reads', weight))
            c_weight = float(read.attributes.get('S.capped', weight))
        else:
            weight, s_weight, e_weight, c_weight = read.weight, read.weight, read.weight, read.weight
        
        if read.strand == 1:
            covrow = covp
            for span in read.junctions():
                l = span[0] - leftmost
                r = span[1] - leftmost
                if l > 0 and r < array_length:
                    block = (l,r)
                    junction_hash = span_to_string(block)
                    J_plus[junction_hash] = J_plus.get(junction_hash, 0) + weight
            
            if read.s_tag:
                pos = read.span[0] - leftmost
                if pos >= 0 and pos < array_length:
                    if read.capped:
                        depth_matrix[Cp, pos] += c_weight
                    else:
                        depth_matrix[Sp, pos] += s_weight
                else:
                    read.s_tag = False
                    read.capped = False
            
            if read.e_tag:
                pos = read.span[1] - leftmost - 1
                if pos >= 0 and pos < array_length:
                    depth_matrix[Ep, pos] += e_weight
                else:
                    read.e_tag = False
        elif read.strand == -1:
            covrow = covm
            for span in read.junctions():
                l = span[0] - leftmost
                r = span[1] - leftmost
                if l > 0 and r < array_length:
                    block = (l,r)
                    junction_hash = span_to_string(block)
                    J_minus[junction_hash] = J_minus.get(junction_hash, 0) + weight
            
            if read.e_tag:
                pos = max(0, read.span[0] - leftmost)
                if pos >= 0 and pos < array_length:
                    depth_matrix[Em, pos] += e_weight
                else:
                    read.e_tag = False
            
            if read.s_tag:
                pos = read.span[1] - leftmost - 1
                if pos >= 0 and pos < array_length:
                    if read.capped:
                        depth_matrix[Cm, pos] += c_weight
                    else:
                        depth_matrix[Sm, pos] += s_weight
                else:
                    read.s_tag = False
        else: # The read has no features other than non-stranded coverage
            covrow = covn
        
        if splice:
            for span in read.ranges:
                l = span[0] - leftmost
                r = span[1] - leftmost
                depth_matrix[covrow, max(l,0):r] += weight
        else:
            l = read.span[0] - leftmost
            r = read.span[1] - leftmost
            depth_matrix[covrow, max(l,0):r] += weight
        
    return depth_matrix, J_plus, J_minus

cpdef str bedgraph(str chrom, int leftmost, np.ndarray depth_matrix, str seqtype='', int strand=0, int sig=1, bint infer_strand=False):
    """Returns a list of bedgraph lines from an array of height values."""
    cdef:
        int p, lastp
        float v, lastv
        str output = ''
        np.ndarray coverage, positions, values, covstranded, strandratio
        bint contiguous
    
    if seqtype.upper() == 'COV' or seqtype == '':
        contiguous = True
        if strand == 1:
            if infer_strand:
                Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
                covstranded = np.sum(depth_matrix[(covp,covm),:],axis=0)
                strandedpositions = np.where(covstranded > 0)[0]
                if strandedpositions.shape[0] > 0: # Some reads to inform the strand
                    strandratio = np.array(np.interp(range(covstranded.shape[0]), strandedpositions, depth_matrix[covp,strandedpositions]/covstranded[strandedpositions]),dtype=np.float32)
                    strandratio[strandratio < .01] = 0
                    strandratio[strandratio > .99] = 1
                    coverage = depth_matrix[covp,:] + depth_matrix[covn,:]*strandratio
                else:
                    strandratio = np.full(depth_matrix.shape[1], .5)
                    coverage = depth_matrix[covn,:]*strandratio
            else:
                coverage = depth_matrix[-3,:]
        elif strand == -1:
            if infer_strand:
                Sp, Ep, Sm, Em, Cp, Cm, covp, covm, covn = range(9)
                covstranded = np.sum(depth_matrix[(covp,covm),:],axis=0)
                strandedpositions = np.where(covstranded > 0)[0]
                if strandedpositions.shape[0] > 0: # Some reads to inform the strand
                    strandratio = np.array(np.interp(range(covstranded.shape[0]), strandedpositions, depth_matrix[covp,strandedpositions]/covstranded[strandedpositions]),dtype=np.float32)
                    strandratio[strandratio < .01] = 0
                    strandratio[strandratio > .99] = 1
                    coverage = depth_matrix[covm,:] + depth_matrix[covn,:]*(1-strandratio)
                else:
                    strandratio = np.full(depth_matrix.shape[1], .5)
                    coverage = depth_matrix[covn,:]*strandratio
            else:
                coverage = depth_matrix[-2,:]
        else:
            coverage = np.sum(depth_matrix[-3:,:], axis=0)
    elif seqtype.upper() == '5P':
        contiguous = False
        if strand == 1:
            coverage = np.sum(depth_matrix[[0,4],:], axis=0)
        elif strand == -1:
            coverage = np.sum(depth_matrix[[2,5],:], axis=0)
        else:
            coverage = np.sum(depth_matrix[[0,4],:], axis=0) - np.sum(depth_matrix[[2,5],:], axis=0)
    elif seqtype.upper() == 'S':
        contiguous = False
        if strand == 1:
            coverage = depth_matrix[0,:]
        elif strand == -1:
            coverage = depth_matrix[2,:]
        else:
            coverage = depth_matrix[0,:] - depth_matrix[2,:]
    elif seqtype.upper() == 'C':
        contiguous = False
        if strand == 1:
            coverage = depth_matrix[4,:]
        elif strand == -1:
            coverage = depth_matrix[5,:]
        else:
            coverage = depth_matrix[4,:] - depth_matrix[5,:]
    elif seqtype.upper() in ['E', '3P']:
        contiguous = False
        if strand == 1:
            coverage = depth_matrix[1,:]
        elif strand == -1:
            coverage = depth_matrix[3,:]
        else:
            coverage = depth_matrix[1,:] - depth_matrix[3,:]
    
    positions = np.where(np.append(coverage[0],np.diff(coverage)) != 0)[0]
    values = coverage[positions]
    lastp, lastv = -1, 0
    for p,v in zip(positions, values):
        if lastp >= 0 and lastv != 0:
            output += '{}\t{}\t{}\t{}\n'.format(chrom, lastp+leftmost, [lastp+1, p][contiguous]+leftmost, round(lastv, sig))
        
        lastp, lastv = p, v
    
    
    if lastv != 0 and lastp >= 0:
        output += '{}\t{}\t{}\t{}\n'.format(chrom, lastp+leftmost, [lastp+1, coverage.shape[0]][contiguous]+leftmost, round(lastv, sig))
    
    return output

cdef str span_to_string((int, int) span):
    """Converts a tuple of two ints to a string connected by ':'"""
    return '{}:{}'.format(span[0], span[1])


cdef (int, int) string_to_span(str string):
    """Converts a string from span_to_string() back into a span"""
    cdef list splitstring = string.split(':')
    return (int(splitstring[0]), int(splitstring[1]))


cdef parse_BED_line(bed_line, chrom_dict, source_dict, source_string=None, s_tag=False, e_tag=False, capped=False, gaps_are_junctions=False, keep_readname=False):
    """Parses one line of a 12- or 15-column BED file into an ELdata namedtuple.
    Examples:
      Ath_chr1 6787 8737 AT1G01020.6 . - 6787 8737 0,0,0 6 282,294,86,90,48,144 0,369,776,1448,1629,1806
      Ath_chr1 6787 8737 . . - 0 0 109,74,116 6 282,294,86,90,48,144 0,369,776,1448,1629,1806 1.0 TAIR10.40 EADADADADADS
    """
    cdef str first, last
    cdef bint condensed
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
        chrom_string, chromStart, end, readname, score, bed_strand = bed_elements[:6]
        blockSizes = str(int(end)-int(chromStart))
        blockStarts = '0'
        label = None
        try:
            weight = float(score)
        except:
            weight = float(1)
    
    if bed_strand == '+':
        strand = 1
    elif bed_strand == '-':
        strand = -1
    else:
        strand = 0
    
    condensed = False
    ranges = get_block_ranges(chromStart, blockStarts, blockSizes)
    if label is not None: # Determine what kind of end labels exist based on the label
        first, last = label[0].upper(), label[-1].upper()
        condensed = first != label[0] or last != label[-1]
        s_tag = e_tag = capped = False
        if strand == 1:
            if first == 'C':
                s_tag = capped = True
            elif first == 'S':
                s_tag = True
            
            if last == 'E':
                e_tag = True
            
            splice = [True if i=='D' else False for i in label[1:-1:2]]
        elif strand == -1:
            if last == 'C':
                s_tag = capped = True
            elif last == 'S':
                s_tag = True
            
            if first == 'E':
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
        source = 0
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight, condensed)

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
    
    flank = genome[chrom][range_start:range_end].upper()
    if strand == -1:
        flank = fu.rc(flank)
    
    return flank

cdef int motif_strand(dict genome, str chrom, int left, int right):
    """pair of splice junction positions based
    on the flanking genomic sequence """
    cdef str flanking_sequence
    cdef set canonical_plus = set(['GTAG','GCAG','ATAC','GAAG'])
    cdef set canonical_minus = set(['CTAC','CTGC','GTAT','CTTC'])
    if chrom not in genome:
        return 0
    
    flanking_sequence = get_flank(genome, chrom, left-1, 1, 'E', 2) + get_flank(genome, chrom, right, 1, 'S', 2)
    flanking_sequence = flanking_sequence.upper()
    if flanking_sequence in canonical_plus:
        return 1
    elif flanking_sequence in canonical_minus:
        return -1
    else:
        return 0

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



cpdef parse_SAM_CIGAR(int pos, list cigartuples, str mdstring, float error_rate=0.1, bint quality_filter=False, int maxindel=10):
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
    cdef bint first_element, jumped, indel_too_large
    cdef str next_match
    ranges = []
    gaps = []
    badgaps = []
    head = 0
    tail = 0
    current_pos = pos
    first_element = True
    jumped = True
    indel_too_large = False
    cigar_len = len(cigartuples)
    mismatches = 0
    exon_size = 0
    i = 0
    mismatch_generator = parse_MD_string(mdstring)
    try:
        match_count = int(next(mismatch_generator))
    except:
        match_count = 0
    
    for i in range(cigar_len):
        operator = cigartuples[i][0]
        o_len = cigartuples[i][1]
        if operator == 4 or operator == 5: # Softclipped ranges are stored in either head or tail
            if first_element:
                head = -1 if operator == 5 else o_len
            else:
                tail = -1 if operator == 5 else o_len
        
        if operator in [0,7,8]: # Match
            if jumped: # The last match was across an intron, dump it
                jumped = False
                if indel_too_large or (quality_filter and exon_size > 0 and mismatches/exon_size > error_rate): # Last exon didn't pass quality threshold
                    del ranges[-1] # Get rid of the last exon which contains the mapping error(s)
                    if len(ranges) == 0: # This removed the leftmost block
                        head = -1 # Hardclip the head
                        gaps = [] # Remove all gaps
                        badgaps = []
                    elif len(ranges) > 0: # Deleted an internal exon; extend the last gap
                        gaps[-1] = (gaps[-2][0], current_pos)
                        badgaps[-1] = True
                        del gaps[-2]
                        del badgaps[-2]
                
                if len(ranges) > 0: # Adjust last exon border to match gap
                    ranges[-1] = (ranges[-1][0], gaps[-1][0])
                    ranges += [(gaps[-1][1],current_pos+o_len)] # Start the next exon from the last junction border
                else:
                    ranges += [(current_pos,current_pos+o_len)] # Add a first exon
                
                mismatches = 0
                exon_size = 0
                indel_too_large = False
            else: # Continuation of exon, update right side
                ranges[-1] = (ranges[-1][0], current_pos+o_len)
            
            exon_size += o_len
            current_pos += o_len
            if operator == 8:
                mismatches += o_len
            else:
                match_count -= o_len
        elif operator == 3: # Skipped region (N)
            leftside = current_pos
            current_pos += o_len
            rightside = current_pos
            if not jumped:
                gaps += [(leftside,rightside)] # Add a new junction
                badgaps += [False]
            else: # No matches encountered in last exon, merge with previous
                gaps[-1] = (gaps[-1][0],rightside)
                badgaps[-1] = True

            jumped = True
        elif operator == 2: # Deletion in query
            current_pos += o_len
            mismatches += o_len
            indel_too_large = indel_too_large or o_len > maxindel
        elif operator == 1: # Insertion in query
            mismatches += o_len
            indel_too_large = indel_too_large or o_len > maxindel
        
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

        first_element = False
    
    if indel_too_large or (quality_filter and exon_size > 0 and mismatches/exon_size > error_rate): # Last exon didn't pass quality threshold
        del ranges[-1]
        if len(gaps) > 0:
            del gaps[-1]
            del badgaps[-1]
        
        tail = -1
    
    assert len(gaps) == len(badgaps)
    return ranges, gaps, badgaps, head, tail


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
        ID_string = string.split(tagsplit)[-1].split('_')[0].upper()
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
    cdef readonly float error_rate
    cdef readonly bint ignore_ends, secondary, remove_noncanonical, remove_gapped_termini, quality_filter
    cdef readonly int max_headclip, sj_shift, max_intron
    def __init__(self, RNAseqDataset dataset, list input_lines):
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
        self.ignore_ends = self.dataset.ignore_ends
        self.secondary = self.dataset.secondary
        self.max_intron = self.dataset.max_intron
        self.remove_noncanonical = self.dataset.remove_noncanonical
        self.remove_gapped_termini = self.dataset.remove_gapped_termini
        self.error_rate = self.dataset.error_rate
        self.quality_filter = self.dataset.quality_filter
        self.max_headclip = self.dataset.max_headclip
        self.sj_shift = self.dataset.sj_shift
    
    cpdef list generate_read(self):
        cdef:
            int s_len, s_tag_len, e_len, e_tag_len, Nmap, counter, input_len, map_number, mate, strand, number_of_blocks
            int i, gap_len, pos, map_strand, junction_strand, chrom_id, start_pos, end_pos, trim_pos, errors, sj_shift, longest_intron
            float weight
            str ID, chrom, seq, aligned_seq, trimmed_nuc
            (bint, bint, int, int) ID_tags = (False, False, 0, 0)
            dict mappings
            list canonical, gaps, ranges, mapping_list, introns, newranges, jstrands, badgaps
            bint stranded, stranded_method, reverse, fiveprime, threeprime, junction_exists, c, b
            (int, int) g
            set intronset
            array.array flankmatch
            RNAseqMapping current_mapping, mate_read
        
        s_tag_len = e_tag_len = 0
        stranded = stranded_method = False
        if self.dataset.stranded: # The read is strand-specific
            stranded_method = True # The method to generate the read is inherently stranded
            reverse = self.dataset.reverse
        else:
            reverse = False
        
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
        
        s_tag_len = max([ID_tags[2], (0, s_len)[self.dataset.s_tag]])
        e_tag_len = max([ID_tags[3], (0, e_len)[self.dataset.e_tag]])
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
            
            mate, strand = self.determine_strand(line, stranded, reverse)
            map_strand = [1, -1][line.is_reverse] * [-1, 1][line.is_read1 or not line.is_paired]
            counter += mate==1
            try:
                map_number = line.get_tag('HI')
            except KeyError:
                map_number = counter
            
            seq = line.query_sequence
            pos = line.reference_start
            chrom_id = line.reference_id
            try:
                chrom = line.header.get_reference_name(chrom_id)
            except:
                chrom = line.reference_name
            
            # Parse the SAM CIGAR string to get mapped positions, splice junction sites, and softclipped positions
            try:
                errors = line.get_tag('nM')
            except KeyError:
                try:
                    errors = line.get_tag('NM')
                except KeyError:
                    errors = 0
            
            try:
                mdstring = line.get_tag('MD')
            except KeyError:
                mdstring = str(len(seq))
            
            ranges, introns, badgaps, head, tail = parse_SAM_CIGAR(pos, line.cigartuples, mdstring, self.error_rate, self.quality_filter)
            number_of_blocks = len(ranges)
            if number_of_blocks == 0: # No exons of passing quality were found
                continue
            
            if tail <= 0:
                aligned_seq = seq[max(0,head):]
            else:
                aligned_seq = seq[max(0,head):-tail]
            
            match_length = len(aligned_seq) - errors
            if match_length < self.dataset.minlen_loose: # Read is short enought to require stringent filtering
                if self.fails_stringent_filters(Nmap, match_length, head, tail, errors):
                    continue
            
            if match_length < 200: # Short read check: Not a homopolymer
                if is_homopolymer(aligned_seq): # Aligned sequence >80% repeat of one nucleotide
                    continue
            
            # EVALUATE SOFTCLIPPED NUCLEOTIDES
            fiveprime = mate == 1
            threeprime = (mate == 1 and not line.is_paired) or mate == 2
            s_tag, e_tag, capped, strand = self.filter_labels_by_softclip(chrom, ranges, seq, s_tag, e_tag, capped, fiveprime, threeprime, strand, head, tail, self.max_headclip)
            if number_of_blocks == 1:
                junction_strand = 0
                canonical = []
            else: # Check if the splice junctions can be used to assign a strand
                junction_strand = self.get_alignment_strand(line, map_strand)
                if strand == 0:
                    strand, canonical = self.evaluate_splice_sites(chrom, junction_strand, introns, badgaps)
                else:
                    strand, canonical = self.evaluate_splice_sites(chrom, strand, introns, badgaps)
                
                if len(self.dataset.genome) > 0: # Use genomic motifs to identify strand
                    jstrands = [motif_strand(self.dataset.genome, chrom, left, right) for left,right in introns]
                    consensus_strand = sum(jstrands)
                    if consensus_strand > 0:
                        strand, consensus_strand = 1, 1
                    elif consensus_strand < 0:
                        strand, consensus_strand = -1, -1
                    
                    canonical = [True if ((js == consensus_strand and consensus_strand !=0) or can) else False for js,can in zip(jstrands, canonical)]

                
                # Check if any splice junctions can be corrected by shifting
                if self.sj_shift > 0:
                    strand, ranges, canonical = self.adjust_splice_sites(strand, ranges, canonical, chrom)
            
            # Reconcile strand information given by start, end, and splice
            junction_exists = sum(canonical) > 0
            if junction_strand != 0 and junction_exists: # At least one strand-informative splice junction exists
                if strand != junction_strand: # Splice disagrees with end tags; remove tags
                    strand = junction_strand
                    s_tag = e_tag = capped = False
            
            if not stranded_method and not s_tag and not e_tag and not junction_exists:
                strand = 0 # No strand information can be found
            
            if not self.remove_noncanonical or (len(self.dataset.genome)==0 and len(self.dataset.sj_set)==0):
                canonical = [True]*len(canonical)
            
            canonical = [c and not b for c,b in zip(canonical, badgaps)] # All bad gaps are False
            if self.remove_gapped_termini:
                ranges, canonical, left_trim, right_trim = self.remove_terminal_gaps(ranges, canonical)
                if strand == 1:
                    s_tag = s_tag and not left_trim
                    capped = capped and not left_trim
                    e_tag = e_tag and not right_trim
                elif strand == -1:
                    s_tag = s_tag and not right_trim
                    capped = capped and not right_trim
                    e_tag = e_tag and not left_trim
            
            # Generate a ReadObject with the parsed attributes above
            read_data = ELdata(chrom_id, 0, strand, ranges, canonical, s_tag, e_tag, capped, round(weight,2), False)
            current_mapping = RNAseqMapping(read_data, attributes = {'errors':errors})
            if current_mapping.get_length() < self.dataset.minlen_strict:
                continue
            
            longest_intron = max([0]+[(j[1]-j[0]) for j in current_mapping.junctions()])
            if longest_intron > self.max_intron:
                # print(f'REMOVING {current_mapping}: intron length {longest_intron}')
                continue

            current_mapping.e_len = e_tag_len
            current_mapping.s_len = s_tag_len
            if map_number not in mappings:
                mappings[map_number] = current_mapping
            else: # merge two mate-pair ReadObjects together
                mate_read = mappings[map_number]
                # if e_tag_added: # Mate alignment cannot extend beyond a trimmed 3' poly(A)
                #     self.update_trimmed_ranges(current_mapping, mate_read, strand)
                
                mate_read.merge(current_mapping)
        
        mapping_list = self.resolve_overlapping_mappings(mappings)
        return mapping_list
    
    cdef (bint, bint, bint, int) filter_labels_by_softclip(self, str chrom, list ranges, str seq, bint s_tag, bint e_tag, bint capped, bint fiveprime, bint threeprime, int strand, int head, int tail, int max_headclip, int eval_length=10):
        """Examines the left and right softclip of the read to determine
        if s_tag, capped, e_tag, or strand attributes should be modified."""
        cdef str leftclip, rightclip
        cdef int left_pos, right_pos, clipstrand
        cdef bint start_plus, start_minus, end_plus, end_minus, capped_plus, capped_minus
        # Get softclipped nucleotide sequences
        if head == 0 and tail == 0: # No softclipping at all, pass defaults
            capped = self.dataset.capped
        else:
            if head < 0:
                s_tag = False if strand == 1 else s_tag
                e_tag = False if strand == -1 else e_tag
                capped = False if strand == 1 else capped
                leftclip = ''
            elif head == 0:
                leftclip = ''
            else:
                leftclip = seq[max(0,head-eval_length):head]
            
            if tail < 0:
                s_tag = False if strand == -1 else s_tag
                e_tag = False if strand == 1 else e_tag
                capped = False if strand == -1 else capped
                rightclip = ''
            elif tail == 0:
                rightclip = ''
            else:
                rightclip = seq[-tail:len(seq)-tail+eval_length]
            
            # Determine strand by clipped sequences
            capped_plus, capped_minus = False, False
            start_plus = head >= 3 and leftclip[-3:] == 'G'*len(leftclip[-3:])
            if start_plus or (strand == 1 and s_tag):
                capped_plus = capped or (0 < len(leftclip) < 5 and leftclip == 'G'*len(leftclip)) or leftclip[-4:] == 'GGGG'
            
            start_minus = tail >= 3 and rightclip[:3] == 'C'*len(rightclip[:3])
            if start_minus or (strand == -1 and s_tag):
                capped_minus = capped or (0 < len(rightclip) < 5 and rightclip == 'C'*len(rightclip)) or rightclip[-4:] == 'CCCC'
            
            end_plus = rightclip.count('A') >= (1.0 - self.dataset.mismatch_rate) * max([1, len(rightclip)])
            end_minus = leftclip.count('T') >= (1.0 - self.dataset.mismatch_rate) * max([1, len(leftclip)])
            
            # Reconcile clipped labels with preset labels
            if (start_plus or end_plus) and not (start_minus or end_minus):
                s_tag = start_plus or (s_tag and 0 <= head < 5)
                e_tag = (end_plus and tail >= 3) or (e_tag and 0 <= tail < 3)
                capped = capped_plus
                strand = 1
            elif (start_minus or end_minus) and not (start_plus or end_plus):
                s_tag = start_minus or (s_tag and 0 <= tail < 5)
                e_tag = (end_minus and head >= 3) or (e_tag and 0 <= head < 3)
                capped = capped_minus
                strand = -1
            else: # softclipping is inconsistent
                if start_plus and end_plus: # both tags agree, override strand
                    strand = 1
                    s_tag, e_tag, capped = start_plus, end_plus, capped_plus
                
                if start_minus and end_minus: # both tags agree, override strand
                    strand = -1
                    s_tag, e_tag, capped = start_minus, end_minus, capped_minus
                
                # default strand wins
                if strand == 1:
                    s_tag = start_plus or (s_tag and 0 <= head < 5)
                    e_tag = end_plus or (e_tag and 0 <= tail < 3)
                    capped = capped_plus
                elif strand == -1:
                    s_tag = start_minus or (s_tag and 0 <= tail < 5)
                    e_tag = end_minus or (e_tag and 0 <= head < 3)
                    capped = capped_minus
        
        # Perform artifact masking against genome sequence
        left_pos = ranges[0][0]
        right_pos = ranges[-1][1]-1
        if self.dataset.has_genome:
            if strand == 1:
                if s_tag:
                    s_tag = not self.matches_masking_sequence(chrom, left_pos, strand, 'S', min([5 if head < 3 else head, len(self.dataset.start_seq)]))
                    capped = capped and s_tag
                if e_tag:
                    e_tag = not self.matches_masking_sequence(chrom, right_pos, strand, 'E', min([eval_length if tail < 3 else tail, len(self.dataset.end_seq)]))
            elif strand == -1:
                if s_tag:
                    s_tag = not self.matches_masking_sequence(chrom, right_pos, strand, 'S', min([5 if tail < 3 else tail,len(self.dataset.start_seq)]))
                    capped = capped and s_tag
                if e_tag:
                    e_tag = not self.matches_masking_sequence(chrom, left_pos, strand, 'E', min([eval_length if tail < 3 else tail, len(self.dataset.end_seq)]))
        
        s_tag = s_tag and fiveprime
        capped = capped and s_tag
        e_tag = e_tag and threeprime
        return s_tag, e_tag, capped, strand
    
    cdef bint matches_masking_sequence(self, str chrom, int position, int strand, str readtype, int length):
        """Evaluates whether a clipped tag matches too closely
        with a genome-templated region could have caused 
        false positive end signal"""
        ## 5'
        flank = get_flank(self.dataset.genome, chrom, position, strand, readtype, length) # Get upstream flanking sequence to start
        
        if len(flank) > 0:
            if readtype == 'S':
                flankmatch = self.dataset.start_array[-length:][::-1]
                flank = flank[::-1]
            elif readtype == 'E':
                flankmatch = self.dataset.end_array[:length]
            
            if fu.oligo_match(fu.nuc_to_int(flank), flankmatch, self.dataset.mismatch_rate):
                return True
        
        return False
    
    def remove_terminal_gaps(self, ranges, canonical, left_trim=False, right_trim=False):
        '''Recursively cuts off first and last exons if they are separated by a
        non-splice gap from the rest of the alignment.'''
        new_left_trim = False
        new_right_trim = False
        if len(canonical) == 0: # No gaps
            return ranges, canonical, left_trim or False, right_trim or False
        elif len(canonical) == 1: # Special case: single gap
            if not canonical[0]: # Pick the longer exon
                biggest = 0
                best_i = 0
                for i in range(len(ranges)):
                    size = ranges[i][1] - ranges[i][0]
                    if size > biggest:
                        best_i = i
                        biggest = size
                
                return [ranges[best_i]], [], left_trim or best_i==1, right_trim or best_i==0
        else:
            if not canonical[0]: # Starting gap
                new_left_trim = True
                ranges = ranges[1:]
                if len(canonical) > 1:
                    canonical = canonical[1:]
                else:
                    canonical = []
            
            if len(canonical) > 0:
                if not canonical[-1]: # Ending gap
                    new_right_trim = True
                    ranges = ranges[:-1]
                    if len(canonical) > 1:
                        canonical = canonical[:-1]
                    else:
                        canonical = []
        
        if new_left_trim or new_right_trim:
            return self.remove_terminal_gaps(ranges, canonical, left_trim or new_left_trim, right_trim or new_right_trim)
        else:
            return ranges, canonical, left_trim, right_trim
    
    cpdef evaluate_splice_sites(self, chrom, int strand, list introns, list badgaps):
        """Returns a list of bools (1 per gap) that summarizes whether
        each splice junction is supported by the available SJDB."""
        cdef list jstrands, canonical
        cdef int l, r, js, junction_strand, consensus_strand
        cdef bint seen_plus, seen_minus
        jstrands = [0]*len(introns)
        seen_plus = False
        seen_minus = False
        for inum in range(len(introns)):
            l,r = introns[inum]
            js = self.check_sjdb(chrom, l, r)
            seen_plus = seen_plus or js == 1
            seen_minus = seen_minus or js == -1
            jstrands[inum] = js
        
        consensus_strand = (strand + sum(jstrands))
        if consensus_strand > 0:
            consensus_strand = 1
        elif consensus_strand < 0:
            consensus_strand = -1
        
        canonical = [True if (js == consensus_strand and consensus_strand !=0) else False for js in jstrands]
        canonical = [c and not b for c,b in zip(canonical, badgaps)] # All bad gaps are False
        # print('evaluate_splice_sites: {}'.format(canonical))
        return consensus_strand, canonical
    
    def adjust_splice_sites(self, strand, ranges, canonical, chrom):
        """Use genome motif information to version of the exon list where borders
        have been shifted to correct poorly-formed splice junctions."""
        cdef int js, i, left, right, shiftleft, shiftright
        # Shift to nearest SJDB first (allow indels)
        if all(canonical):
            return strand, ranges, canonical
        
        if self.dataset.has_genome or self.dataset.sj_set is not None:
            for i in range(1, len(ranges)): # Iterate over introns
                if not canonical[i-1]: # Skip any that are canonical
                    left = ranges[i-1][1]
                    right = ranges[i][0]
                    shiftleft, shiftright, js = self.shift_junction(chrom, left, right, strand)
                    if strand == 0 or js == strand:
                        strand = js
                        ranges[i-1] = (ranges[i-1][0], shiftleft)
                        ranges[i] = (shiftright, ranges[i][1])
                        canonical[i-1] = True
        
        return strand, ranges, canonical
    
    cdef (int, int, int) shift_junction(self, str chrom, int left, int right, int strand, bint indel=True):
        """Returns a new left/right border of an intron and its strand."""
        cdef int i, js
        for i in range(0, self.sj_shift+1): # Check SJDB, then genome
            if strand >= 0:
                if (chrom, left+i, right+i, 1) in self.dataset.sj_set:
                    # print("Found forward SJDB at shift +{}".format(i))
                    return (left+i, right+i, 1)
                elif (chrom, left-i, right-i, 1) in self.dataset.sj_set:
                    # print("Found forward SJDB at shift -{}".format(i))
                    return (left-i, right-i, 1)
            
            if strand <= 0:
                if (chrom, left+i, right+i, -1) in self.dataset.sj_set:
                    # print("Found reverse SJDB at shift +{}".format(i))
                    return (left+i, right+i, -1)
                elif (chrom, left-i, right-i, -1) in self.dataset.sj_set:
                    # print("Found reverse SJDB at shift -{}".format(i))
                    return (left-i, right-i, -1)
            
            js = motif_strand(self.dataset.genome, chrom, left+i, right+i)
            if (js == strand and strand != 0) or (js !=0 and strand == 0):
                # print("Found {} motif at shift +{}".format(js, i))
                return (left+i, right+i, js)
            
            js = motif_strand(self.dataset.genome, chrom, left-i, right-i)
            if (js == strand and strand != 0) or (js !=0 and strand == 0):
                # print("Found {} motif at shift -{}".format(js, i))
                return (left-i, right-i, js)
    
    cdef list resolve_overlapping_mappings(self, dict mappings):
        """Decision tree for simplifying multimappers """
        cdef RNAseqMapping mapping, v, previous_mapping
        cdef list output_mappings, ordered_mappings
        cdef int length
        if len(mappings.values()) <= 1:
            output_mappings = list(mappings.values())
        else:
            output_mappings = []
            ordered_mappings = [mapping for length,mapping in sorted([(v.span[1]-v.span[0],v) for v in mappings.values()])]
            for mapping in ordered_mappings:
                for previous_mapping in output_mappings:
                    if mapping.overlaps(previous_mapping) and previous_mapping.attributes.get('errors',0) <= mapping.attributes.get('errors', 0):
                        continue
                    
                output_mappings.append(mapping)
        
        return output_mappings
    
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
    
    cdef (int, int) determine_strand(self, line, bint stranded, bint reverse=False):
        """Determine the RNA strand of a pysam object"""
        cdef int mate, strand
        mate = int(not(line.is_read1 or line.is_read2)) + line.is_read1 + 2*line.is_read2
        strand = stranded * (((mate==1) ^ (reverse ^ line.is_reverse))*2 - 1)
        return mate, strand
    
    cdef bint should_skip(self, line):
        """The read should not be processed."""
        return line.is_unmapped or ((line.is_supplementary or line.is_secondary) and not self.secondary) or (line.is_paired and not line.is_proper_pair)
    
    cdef bint fails_stringent_filters(self, int Nmap, int match_length, int head, int tail, int errors):
        """Reads below the 'minlen_loose' length should be treated
        more stringently: no allowed softclipping, multimapping, or mismatches.
        Absolutely require the length to be longer than minlen_strict."""
        return (head > 0 or tail > 0 or errors > 0 or Nmap > 1) and match_length < self.dataset.minlen_strict
    
    cdef int get_alignment_strand(self, line, int strand):
        """Returns 1(+), -1(-), or 0(.) if one of the BAM
        splice tags (XS, ts) contains strand information."""
        cdef str js
        cdef int alignment_strand
        alignment_strand = strand
        try:
            js = line.get_tag('XS')
            if js == '+':
                alignment_strand = 1
            elif js == '-':
                alignment_strand = -1
        except KeyError:
            try:
                js = line.get_tag('ts')
                alignment_strand = {'+':strand, '-':-strand}[js]
            except KeyError:
                alignment_strand = 0
        
        return alignment_strand
    
    cdef int check_sjdb(self, str chrom, int start, int end):
        '''Returns the strand of'''
        cdef int strand
        strand = 0
        strand += int((chrom, start, end, 1) in self.dataset.sj_set)
        strand -= int((chrom, start, end, -1) in self.dataset.sj_set)
        return strand


def read_generator(fileconn, RNAseqDataset dataset, str file_type, int max_gap, float minimum_proportion, int max_intron=100000):
    """Yields a contiguous chunk of reads from the input file
    separated on either side by a gaps > max_gap"""
    cdef RNAseqMapping read, outread
    cdef int l, r, old_chrom, old_l, old_r, rightmost, k
    cdef (int, int) j
    cdef float read_weight, span_weight, current_cov
    cdef set covered_positions
    cdef list passed_positions
    if file_type in ['elr','elr.gz']:
        add_read = dataset.add_read_from_ELR
    elif file_type == 'bed':
        add_read = dataset.add_read_from_BED
    elif file_type == 'bam' or file_type == 'sam':
        add_read = dataset.add_read_from_BAM
    else:
        return
    
    end_positions = Counter() # Keep track of where reads end to maintain a tally of coverage depth
    old_chrom, old_l, old_r, rightmost, span_start = -1, -1, -1, -1, -1
    current_cov = 0
    span_weight = 0
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
        if len(dataset.read_list) == 0:
            continue
        
        read = dataset.read_list[-1]
        if read.is_malformed() or any([(j[1]-j[0]) > max_intron for j in read.junctions()]):
            del dataset.read_list[-1]
            continue
        
        l, r = read.span
        read_weight = read.weight * (r-l)
        current_cov += read.weight
        if old_chrom == -1: # Uninitialized; add the read and make no other decisions
            span_start = l
            span_length = r - span_start
            rightmost = r
        elif read.chrom != old_chrom or l >= rightmost + max_gap: # The last locus is definitely finished; dump the read list
            if read.chrom != old_chrom:
                yield [outread for outread in dataset.read_list[:-1]]
            else:
                yield [outread for outread in dataset.read_list[:-1] if outread.span[1] < l]
            
            dataset.read_list = [read]
            span_start = l
            span_weight = 0
            current_cov = read.weight
            end_positions = Counter()
            rightmost = r
        elif l > old_l: # Read advanced, but not by enough to automatically cut
            passed_positions = [k for k in end_positions.keys() if k <= l]
            for k in passed_positions:
                current_cov -= end_positions.pop(k)
            
            if current_cov * span_length < minimum_proportion * span_weight: # Current cov is sufficiently lower than mean cov to cause a break
                yield [outread for outread in dataset.read_list[:-1] if outread.span[1] < l]
                dataset.read_list = [read]
                span_start = l
                span_weight = 0
                current_cov = read.weight
                end_positions = Counter()
                rightmost = r
        
        end_positions[r] += read.weight # Add the read's weight to the position where the read ends
        span_weight += read_weight
        if r > rightmost: rightmost = r
        span_length = rightmost - span_start
        old_chrom, old_l, old_r = read.chrom, l, r
    
    # Dump the remaining reads
    yield dataset.read_list
    fileconn.close()

def generate_subchunks(list list_of_reads, list split_positions):
    cdef:
        RNAseqMapping read
        int lasti, i, sp, r
        set ignore
    
    lasti = 0
    position = iter(split_positions)
    sp = next(position)
    ignore = set()
    for i,read in enumerate(list_of_reads):
        if read.span[0] > sp: # The current read passes the split
            if i > lasti:
                yield [list_of_reads[r] for r in range(lasti, i) if r not in ignore]
            
            lasti = i
            try:
                while sp < read.span[0]:
                    sp = next(position)
            except StopIteration:
                lasti = i
                break
        
        if read.span[1] > sp:
            ignore.add(i)
    
    yield list_of_reads[lasti:]


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

cpdef bint has_ends(list list_of_reads, bint require_cap):
    cdef RNAseqMapping read
    cdef bint sp, cp, ep, sm, cm, em
    sp, cp, ep, sm, cm, em = False, False, False, False, False, False
    for read in list_of_reads:
        if read.strand == 1:
            sp = sp or read.s_tag
            cp = cp or read.capped
            ep = ep or read.e_tag
            if (not require_cap or cp) and sp and ep:
                return True
        elif read.strand == -1:
            sm = sm or read.s_tag
            cm = cm or read.capped
            em = em or read.e_tag
            if (not require_cap or cm) and sm and em:
                return True
        
    return False

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
                gap_left = i+1
            
            gap_length += 1
            gap_is_long_enough = gap_length > maxgap
        else: # Not in a gap
            if gap_is_long_enough: # The last position was the end of a gap
                gap_right = i
                gaps.append((gap_left, gap_right))
            
            gap_length = 0
            gap_is_long_enough = False
    
    return gaps

