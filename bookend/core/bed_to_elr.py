#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import bookend.core.cython_utils._rnaseq_utils as ru
if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import bed_to_elr_parser as parser

class BEDtoELRconverter:
    def __init__(self, args):
        """Converts each line of BED-formatted input to ELR"""
        self.input = args['INPUT']
        self.chrom_file = args['CHROMS']
        self.output = args['OUTPUT']
        self.header = args['HEADER']
        self.source = args['SOURCE']
        self.junctions = args['JUNCTIONS']
        self.start = args['START']
        self.capped = args['CAPPED']
        self.end = args['END']
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            self.output_file = open(self.output, 'w')
        
        if self.header is None:
            self.header_file = self.output_file
        else:
            self.header_file = open(self.header, 'w')
        
        if self.source:
            source_array = [self.source]
        else:
            source_array = None
        
        self.chrom_array = None
        if self.chrom_file:
            self.chrom_array = [l.strip() for l in open(self.chrom_file, 'r')]
        
        self.linecount = 0
        self.readcount = 0
        self.dataset = ru.RNAseqDataset(source_array=source_array, chrom_array=self.chrom_array)
    
    def process_input(self):
        """Yield ELR lines from each line of a BED input."""
        current_source_index = self.dataset.source_index
        current_chrom_index = self.dataset.chrom_index
        
        if current_chrom_index > 0: # Header already exists
            for i,c in enumerate(self.dataset.chrom_array):
                self.output_line('#C {} {}'.format(i,c), self.header_file)
        
        if current_source_index > 0:
            for i,s in enumerate(self.dataset.source_array):
                self.output_line('#S {} {}'.format(i,s), self.header_file)
        
        
        bed_in = open(self.input, 'r')

        for bed_line in bed_in:
            if bed_line[0] == '#': # A header is being passed from the BED file
                header_line = bed_line.rstrip().split(' ')
                if header_line[0] == '#S':
                    self.dataset.add_source(header_line[-1])
                if header_line[0] == '#C':
                    self.dataset.add_chrom(header_line[-1])
                
                self.output_line(bed_line.rstrip(), self.header_file)
                continue
            
            current_chrom_index = self.dataset.chrom_index
            current_source_index = self.dataset.source_index
            self.dataset.add_read_from_BED(bed_line, source_string=self.source, s_tag=self.start, e_tag=self.end, capped=self.capped, gaps_are_junctions=self.junctions)
            self.linecount += 1
            self.readcount += self.dataset.read_list[-1].weight
            if self.dataset.chrom_index > current_chrom_index: # A new chromosome was encountered
                self.output_line('#C {} {}'.format(len(self.dataset.chrom_array)-1, self.dataset.chrom_array[-1]), self.header_file)
                current_chrom_index = self.dataset.chrom_index
            
            if self.dataset.source_index > current_source_index: # A new source was encountered
                self.output_line('#S {} {}'.format(len(self.dataset.source_array)-1, self.dataset.source_array[-1]), self.header_file)
                current_source_index = self.dataset.source_index
            
            while len(self.dataset.read_list) > 0:
                out_line = self.dataset.pop_read()
                self.output_line(out_line, self.output_file)
        
        if self.output != 'stdout':
            self.output_file.close()
        
        if self.header_file != self.output_file and self.header_file != 'stdout':
            self.header_file.close()

    def output_line(self, line, output_file):
        """Takes a list of bed lines and writes
        them to the output stream.
        """
        if output_file == 'stdout':
            print(line)
        else:
            output_file.write('{}\n'.format(line.rstrip()))
    
    def run(self):
        print(self.display_options())
        if not self.input.split('.')[-1].lower() in ['bed','bed12']:
            print("ERROR: input must be in the BED12 format (.bed, .bed12)")
            return 1
        
        self.process_input()
        print(self.display_summary())
        return 0
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend bed-to-elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        options_string += "  Header output: {}\n".format(self.header)
        options_string += "  Output source: {}\n".format(self.source)
        options_string += "  *** Parameters ***\n"
        options_string += "  All gaps are junctions (-j): {}\n".format(self.junctions)
        options_string += "  All 5' ends are starts (-s): {}\n".format(self.start)
        options_string += "  All 5' ends are capped (-c): {}\n".format(self.capped)
        options_string += "  All 3' ends are polyA (-e):  {}\n".format(self.end)
        return options_string
    
    def display_summary(self):
        summary_string = ''
        summary_string += 'Processed {} lines ({} total reads).\n'.format(self.linecount, round(self.readcount,2))
        return summary_string
    

if __name__ == '__main__':
    args = vars(parser.parse_args())
    obj = BEDtoELRconverter(args)
    sys.exit(obj.run())
        
