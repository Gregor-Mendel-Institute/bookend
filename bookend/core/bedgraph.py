#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import time
import gzip
from pysam import AlignmentFile
import bookend.core.cython_utils._rnaseq_utils as ru
from bookend.core.elr_combine import ELRcombiner

if __name__ == '__main__':
    sys.path.append('../../bookend')


class Bedgrapher:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.start_time = time.time()
        self.type = args['TYPE'].upper()
        self.strand = args['STRAND']
        self.intstrand = {'.':0,'+':1,'-':-1}[self.strand]
        self.output = args['OUT']
        self.input = args['INPUT']    
        # self.scale = args['SCALE']    
        if len(self.input) == 1:
            self.input = self.input[0]
            if self.input_is_valid(self.input):
                self.file_type = self.file_extension(self.input)
                if self.file_type in ['bam','sam']:
                    try:
                        self.input_file = AlignmentFile(self.input)
                    except:
                        print("ERROR: BAM/SAM input must be position-sorted and indexed.")
                        sys.exit(1)
                    
                    self.dataset = ru.RNAseqDataset(chrom_array=self.input_file.header.references)
                elif self.file_type == 'elr.gz':
                    self.dataset = ru.RNAseqDataset()
                    self.input_file = gzip.open(self.input, 'rt')
                else:
                    self.dataset = ru.RNAseqDataset()
                    self.input_file = open(self.input,'r')
            else:
                print("\nERROR: input file must be a valid format (BED, ELR, BAM, SAM).")
                sys.exit(1)
        elif len(self.input) > 1: # Interleave multiple input files for assembly
            if not all([self.input_is_valid(filename, valid_formats=['elr','elr.gz']) for filename in self.input]):
                print("\nERROR: Multi-input assembly can only be performed on position-sorted ELR files.")
                sys.exit(1)
            
            self.dataset = ru.RNAseqDataset()
            self.file_type = 'elr'
            combine_args = {
                'INPUT':self.input,
                'OUTPUT':'stdout',
                'TEMPDIR':'{}_combinetmp'.format(self.input[0])
            }
            combiner = ELRcombiner(combine_args)
            self.input_file = combiner.combine_files(combiner.input, combiner.output_file, iterator=True)
        else:
            print("\nERROR: No input file(s) provided.")
            sys.exit(1)
        
        self.generator = ru.read_generator(self.input_file, self.dataset, self.file_type, 0, 0)
        self.chunk_counter = 0
        self.output_file = open(self.output,'w')
    
    def process_entry(self, chunk):
        if len(chunk) > 0:
            chrom = self.dataset.chrom_array[chunk[0].chrom]
            self.chunk_counter += 1
            leftmost, rightmost = ru.range_of_reads(chunk)
            depth_matrix, J_plus, J_minus = ru.build_depth_matrix(leftmost, rightmost, tuple(chunk), use_attributes=True, splice=True)
            self.output_file.write(ru.bedgraph(chrom, leftmost, depth_matrix, self.type, self.intstrand))
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend bedgraph |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                                       {}\n".format(self.input)
        options_string += "  Output file (-o):                                 {}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Bedgraph read type filter (--type):               {}\n".format(self.type)
        options_string += "  Bedgraph output strand filter (--strand):         {}\n".format(self.strand)
        return options_string
    
    def display_summary(self):
        summary = ''
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
    
    def filter_readtypes(self, chunk):
        """Given the output label type, filter the chunk for only those reads."""
        if self.type in ['','COV']:
            return chunk
        elif self.type == '5P':
            return [read for read in chunk if read.s_tag]
        elif self.type == 'S':
            return [read for read in chunk if read.s_tag and not read.capped]
        elif self.type == 'C':
            return [read for read in chunk if read.capped]
        elif self.type in ['3P','E']:
            return [read for read in chunk if read.e_tag]
    
    def run(self):
        """Executes end labeling on all reads."""
        if self.output != 'stdout':
            print(self.display_options())
        
        for chunk in self.generator:
            chunk = self.filter_readtypes(chunk)
            self.process_entry(chunk)
        
        if len(self.input) == 1:
            self.output_file.close()
        
        self.end_time = time.time()
        print(self.display_summary())

if __name__ == '__main__':
    from argument_parsers import assemble_parser as parser
    args = vars(parser.parse_args())
    obj = Bedgrapher(args)
    sys.exit(obj.run())

