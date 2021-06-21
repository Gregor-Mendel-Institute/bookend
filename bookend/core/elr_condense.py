#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import time
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._assembly_utils as au
from bookend.core.elr_sort import ELRcombiner

if __name__ == '__main__':
    sys.path.append('../../bookend')

class Condenser:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.start_time = time.time()
        self.output = args['OUT']
        self.max_gap = args['MAX_GAP']
        self.end_cluster = args['END_CLUSTER']
        self.min_overhang = args['MIN_OVERHANG']
        self.min_cov = args['MIN_COV']
        self.min_intron_length = args['MIN_INTRON_LEN']
        self.intron_filter = args['INTRON_FILTER']
        self.min_proportion = args['MIN_PROPORTION']
        self.cap_bonus = args['CAP_BONUS']
        self.minlen = args['MINLEN']
        self.input = args['INPUT']
        self.antisense_filter = 0.001
        if self.input_is_valid(self.input):
            self.file_type = self.file_extension(self.input)
            self.dataset = ru.RNAseqDataset()
            self.input_file = open(self.input)
        else:
            print("\nERROR: input file must be a valid format (BED, ELR, BAM, SAM).")
            sys.exit(1)
        
        self.output_type = 'elr'
        self.generator = ru.read_generator(self.input_file, self.dataset, self.file_type, self.max_gap, 0)
        self.chunk_counter = 0
        self.transcripts_written = 0
        self.bases_used = 0
        self.output_temp = open(self.output+'.tmp', 'w')
        #self.output_file = open(self.output,'w')
    
    def sort_output(self):
        self.output_temp.close()
        sorter = ELRsorter({'INPUT':self.output+'.tmp','OUT':self.output,'FORCE':True})
        sorter.run()
    
    def output_transcripts(self, transcript):
        """Writes the RNAseqMapping object 'transcript' to an output stream,
        formatted as output_type."""
        output_line = transcript.write_as_elr()
        self.output_temp.write(output_line+'\n')

    def process_entry(self, chunk):
        if len(chunk) == 1:
            self.chunk_counter += 1
            transcript = chunk[0]
            if self.passes_all_checks(transcript):
                self.output_transcripts(transcript)
                self.transcripts_written += 1
                self.bases_used += transcript.weight*transcript.get_length()
        elif len(chunk) > 1:
            chrom = chunk[0].chrom
            self.chunk_counter += 1
            locus = au.Locus(
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
                intron_filter=self.intron_filter, 
                ignore_ends=False,
                allow_incomplete=True,
                require_cap=False,
                simplify=False
            )
            self.chunk_counter = locus.chunk_number
            total_bases = locus.bases
            if total_bases > 0:
                for transcript in locus.transcripts:
                    if self.passes_all_checks(transcript):
                        self.output_transcripts(transcript)
                        self.transcripts_written += 1
                        self.bases_used += transcript.attributes['bases']
            
            del locus
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend assemble |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                                       {}\n".format(self.input)
        options_string += "  Output file (-o):                                 {}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Max allowed gap in coverage (--max_gap):          {}\n".format(self.max_gap)
        options_string += "  Max end cluster distance (--end_cluster):         {}\n".format(self.end_cluster)
        options_string += "  Min spanning bases (--min_overhang):              {}\n".format(self.min_overhang)
        options_string += "  *** Filters ***\n"
        options_string += "  Min bp transcript length (--min_len):              {}\n".format(self.minlen)
        options_string += "  Min isoform coverage (--min_cov):                 {}\n".format(self.min_cov)
        options_string += "  Min isoform contribution (--min_proportion):      {}\n".format(self.min_proportion)
        options_string += "  Min retained intron proportion (--intron_filter): {}\n".format(self.intron_filter)
        return options_string
    
    def display_summary(self):
        summary = 'Wrote {} transcripts ({} bases).'.format(self.transcripts_written, self.bases_used)
        summary += 'Total elapsed time: {}'.format(round(self.end_time - self.start_time, 5))
        return summary
    
    def file_extension(self, filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        else:
            return split_name[-1].lower()
    
    def input_is_valid(self, filename, valid_formats=['elr']):
        """Boolean if the file is a format that Assembler can parse."""
        if self.file_extension(filename) in valid_formats:
            return True
        else:
            return False
    
    def passes_all_checks(self, transcript):
        """Determines if a transcript model passes the collection
        of filters set by the commandline arguments.
        'transcript' is an RNAseqMapping object,
        see _rnaseq_utils.pyx for details.
        """
        if transcript.coverage < self.min_cov: return False
        if transcript.attributes['length'] < self.minlen: return False
        return True
    
    def run(self):
        """Executes end labeling on all reads."""
        print(self.display_options())
        wrote_header = self.output_type == 'gtf'
        for locus in self.generator:
            if not wrote_header:
                self.output_temp.write('\n'.join(self.dataset.dump_header())+'\n')
                wrote_header = True
            
            self.process_entry(locus)
        
        self.sort_output()
        self.end_time = time.time()
        print(self.display_summary())

if __name__ == '__main__':
    from argument_parsers import condense_parser as parser
    args = vars(parser.parse_args())
    obj = Condenser(args)
    sys.exit(obj.run())

