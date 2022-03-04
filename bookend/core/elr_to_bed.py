#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import gzip
from bookend.core.cython_utils._rnaseq_utils import RNAseqDataset
if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import bed_to_elr_parser as parser

class ELRtoBEDconverter:
    def __init__(self, args):
        """Converts each line of BED-formatted input to ELR"""
        self.input = args['INPUT']
        self.output = args['OUTPUT']
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            self.output_file = open(self.output, 'w')
        
        self.linecount = 0
        self.readcount = 0
        self.dataset = RNAseqDataset()

    def process_input(self):
        """Yield a BED line from each line of a ELR input."""
        if self.input.lower().endswith('.elr.gz'):
            elr_in = gzip.open(self.input, 'rt')
        else:
            elr_in = open(self.input, 'r')
        
        for elr_line in elr_in:
            if elr_line[0] == '#':
                header_line = elr_line.rstrip().split(' ')
                if header_line[0] == '#S':
                    self.dataset.add_source(header_line[-1])
                if header_line[0] == '#C':
                    self.dataset.add_chrom(header_line[-1])
                
                continue

            self.dataset.add_read_from_ELR(elr_line)
            self.linecount += 1
            while len(self.dataset.read_list) > 0:
                self.readcount += self.dataset.read_list[-1].weight
                bed_line = self.dataset.pop_read('bed')
                self.output_line(bed_line, self.output_file)
            
            current_source_index = self.dataset.source_index
            current_chrom_index = self.dataset.chrom_index
        
        if self.output != 'stdout':
            self.output_file.close()
    
    def output_line(self, line, output_file):
        """Takes a list of bed lines and writes
        them to the output stream.
        """
        if output_file == 'stdout':
            print(line)
        else:
            output_file.write('{}\n'.format(line.rstrip()))
    
    def run(self):
        if __name__ == '__main__':print(self.display_options())
        if not (self.input.lower().endswith('.elr') or self.input.lower().endswith('.elr.gz')):
            print("ERROR: input must be in the ELR format (.elr)")
            return 1
        
        self.process_input()
        print(self.display_summary())
        return 0
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend elr-to-bed |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        options_string += "  Header output: {}\n".format(self.header)
        return options_string
    
    def display_summary(self):
        summary_string = ''
        summary_string += 'Processed {} lines ({} total reads).\n'.format(self.linecount, round(self.readcount,2))
        return summary_string
    

if __name__ == '__main__':
    args = vars(parser.parse_args())
    obj = ELRtoBEDconverter(args)
    sys.exit(obj.run())
