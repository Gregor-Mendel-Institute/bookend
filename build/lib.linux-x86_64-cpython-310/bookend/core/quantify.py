#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import time
import gzip
import pysam
import copy
from multiprocessing.dummy import Pool as ThreadPool
import bookend.core.cython_utils._rnaseq_utils as ru
from bookend.core.cython_utils._assembly_utils import Locus
from bookend.core.elr_combine import ELRcombiner

if __name__ == '__main__':
    sys.path.append('../../bookend')


class Quantifier:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        print(args)
        self.start_time = time.time()
        self.input = args['INPUT']
        self.output = args['OUT']
        self.reference = args['REFERENCE']
        self.max_gap = args['MAX_GAP']
        self.end_cluster = args['END_CLUSTER']
        self.max_intron = args['MAX_INTRON']
        self.min_overhang = args['MIN_OVERHANG']
        self.verbose = args['VERBOSE']
        self.antisense_filter = 0.001

        self.gene_attr = args['PARENT_ATTR_GENE']
        self.gene_attr_child = args['CHILD_ATTR_GENE']
        self.gtf_parent = args['GFF_PARENT']
        self.gtf_child = args['GFF_CHILD']
        self.gff_parent = args['GFF_PARENT']
        self.gff_child = args['GFF_CHILD']
        self.refid_parent = args['PARENT_ATTR_TRANSCRIPT']
        self.refid_child = args['CHILD_ATTR_TRANSCRIPT']
        self.gene_delim = args['GENE_DELIM']

        if self.file_extension(self.output) != 'tsv':
            self.output += '.tsv'
        
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        config_defaults['reference'] = self.reference
        self.dataset = ru.RNAseqDataset(
            genome_fasta=None, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        
        if len(self.input) == 1:
            self.input = self.input[0]
            if self.input_is_valid(self.input, valid_formats=['elr','elr.gz']):
                self.file_type = self.file_extension(self.input)
                if self.file_type == 'elr.gz':
                    self.input_file = gzip.open(self.input, 'rt')
                else:
                    self.input_file = open(self.input, 'r')
            else:
                print("\nERROR: input file must be a valid ELR file.")
                sys.exit(1)
        elif len(self.input) > 1: # Interleave multiple input files for assembly
            if not all([self.input_is_valid(filename, valid_formats=['elr','elr.gz']) for filename in self.input]):
                print("\nERROR: Multi-input assembly can only be performed on position-sorted ELR files.")
                sys.exit(1)
            
            self.file_type = 'elr'
            combine_args = {
                'INPUT':self.input,
                'OUTPUT':'stdout',
                'TEMPDIR':'{}_combinetmp'.format(self.input[0])
            }
            combiner = ELRcombiner(combine_args)
            self.input_file = combiner.combine_files(combiner.input, combiner.output_file, iterator=True)
            header = next(self.input_file)
            self.dataset.chrom_array = combiner.dataset.chrom_array
            self.dataset.source_array = combiner.dataset.source_array
            self.dataset.chrom_dict = combiner.dataset.chrom_dict
            self.dataset.source_dict = combiner.dataset.source_dict
        else:
            print("\nERROR: No input file(s) provided.")
            sys.exit(1)
        
        self.chunk_counter = 0
        self.transcript_counter = 0
        self.total_bases = 0
    
    def generator(self):
        '''Returns an iterator that leafs together
        reads and reference transcripts.'''
        self.read_generator = ru.read_generator(self.input_file, self.dataset, self.file_type, self.max_gap, 0, self.max_intron)
        readchunk = next(self.read_generator)
        for chrom in self.dataset.chrom_array:
            self.transcript_generator = ru.read_iterator(self.dataset.reference_dict[chrom], self.dataset, self.max_gap, minimum_proportion=0, max_intron=10^9)
            refchunk = next(self.transcript_generator)
        

        

    def make_config_dicts(self):
        """Converts commandline input into three config dicts
        to pass to the RNAseqDataset."""
        config_defaults = copy.copy(ru.config_defaults)
        gtf_defaults = copy.copy(ru.gtf_defaults)
        gff_defaults = copy.copy(ru.gff_defaults)
        config_defaults['min_reps'] = 0
        config_defaults['cap_percent'] = 0
        config_defaults['verbose'] = self.verbose
        config_defaults['gene_delim'] = self.gene_delim
        if self.gtf_parent: gtf_defaults['parent_types'] = set(self.gtf_parent)
        if self.gtf_child: gtf_defaults['child_types'] = set(self.gtf_child)
        if self.gff_parent: gff_defaults['parent_types'] = set(self.gff_parent)
        if self.gff_child: gff_defaults['child_types'] = set(self.gff_child)
        if self.gene_attr is None:
            self.gene_attr = gtf_defaults['parent_key_gene']
        else:
            gff_defaults['parent_key_gene'] = self.gene_attr
            gtf_defaults['parent_key_gene'] = self.gene_attr
        
        if self.gene_attr_child is not None:
            gff_defaults['child_key_gene'] = self.gene_attr_child
        
        if self.refid_parent is not None:
            gff_defaults['parent_key_transcript'] += self.refid_parent
            gtf_defaults['parent_key_transcript'] += self.refid_parent
        
        if self.refid_child is not None:
            gff_defaults['child_key_transcript'] += self.refid_child
            gtf_defaults['child_key_transcript'] += self.refid_child
        
        return config_defaults, gtf_defaults, gff_defaults

    def assign_reads_to_reference(self, reads, reference):
        '''Given a list of reads and a list of reference transcripts,
        create a locus with transcripts as full paths, and assign
        read weights to each path. Return an ordered list of bases assigned
        to each transcript and unassigned.'''
        pass


    def process_entry(self, chunk):
        STOP_AT=float('inf')
        # STOP_AT=1000000
        if len(chunk) > 0:            
            chrom = chunk[0].chrom
            self.chunk_counter += 1
            if self.verbose:
                print("\nProcessing chunk {} ({}:{}-{}, {} reads)".format(
                    self.chunk_counter,
                    self.dataset.chrom_array[chrom], chunk[0].left(),chunk[-1].right(),
                    len(chunk)
                ), end=" ")
            
            locus = Locus(
                chrom=chrom, 
                chunk_number=self.chunk_counter, 
                list_of_reads=chunk,
                max_gap=self.max_gap,
                end_cluster=self.end_cluster,
                min_overhang=self.min_overhang, 
                reduce=True, 
                minimum_proportion=self.min_proportion, 
                min_intron_length=self.min_intron_length, 
                antisense_filter=self.antisense_filter, 
                cap_bonus=self.cap_bonus,
                cap_filter=self.cap_filter, 
                complete=False, 
                verbose=self.verbose, 
                naive=self.ignore_sources, 
                intron_filter=self.intron_filter,
                splittable=True,
                ignore_ends=self.ignore_labels, 
                allow_incomplete=self.incomplete,
                require_cap=self.require_cap,
                min_start=self.min_start,
                min_end=self.min_end,
                assemble=False,
                truncation_filter=self.truncation_filter,
            )
            self.chunk_counter = locus.chunk_number
            total_bases = locus.bases
            transcripts_written = 0
            self.covfile.write('{}\t{}\n'.format(transcript.attributes['transcript_id'], '\t'.join([str(round(v,1)) for v in source_cov])))
            if self.verbose:
                print('{} transcripts from {}/{} bases ({}%)'.format(
                    transcripts_written, round(bases_used,1), round(total_bases,1), round(bases_used/total_bases*100,2)), end=" ")
            
            if chunk[0].left() >= STOP_AT:
                sys.exit()
            
            self.transcript_counter += transcripts_written
            
            del locus
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend assemble |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input files:                                      {}\n".format(self.input)
        options_string += "  Output file (-o):                                 {}\n".format(self.output)
        options_string += "  Reference file (-o):                              {}\n".format(self.reference)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Max allowed gap in coverage (--max_gap):          {}\n".format(self.max_gap)
        options_string += "  Max end cluster distance (--end_cluster):         {}\n".format(self.end_cluster)
        options_string += "  Min spanning bases (--min_overhang):              {}\n".format(self.min_overhang)
        return options_string
    
    def display_summary(self):
        summary = ''
        summary += '\nTotal elapsed time: {}'.format(round(self.end_time - self.start_time, 5))
        return summary
    
    def file_extension(self, filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        elif split_name[-1].lower() == 'gz':
            return '.'.join([n.lower() for n in split_name[-2:]])
        else:
            return split_name[-1].lower()
    
    def input_is_valid(self, filename, valid_formats=['bed','elr','bam','sam','gtf','gff3','gff','elr.gz']):
        """Boolean if the file is a format that Assembler can parse."""
        if self.file_extension(filename) in valid_formats:
            return True
        else:
            return False
    
    def run(self):
        """Executes end labeling on all reads."""
        if self.output == 'stdout':
            self.output_file = sys.stdout
        else:
            print(self.display_options())
            self.output_file=open(self.output, 'w')

        self.output_file.write('{}\n'.format('\t'.join(self.dataset.source_array)))
        for chunk in self.generator:
            self.process_entry(chunk)
        
        if len(self.input) == 1:
            self.output_file.close()

        self.end_time = time.time()
        print(self.display_summary())

if __name__ == '__main__':
    from argument_parsers import quantify_parser as parser
    args = vars(parser.parse_args())
    obj = Quantifier(args)
    sys.exit(obj.run())

