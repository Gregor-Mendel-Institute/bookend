#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

class ELRsorter:
    def __init__(self, args):
        """Converts each line of BED-formatted input to ELR"""
        self.input = args['INPUT']
        self.output = args['OUT']
        self.force = args['FORCE']
        self.read_tuples = []
        self.linecount = 0
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            if os.path.exists(self.output) and not self.force:
                print("ERROR: output file already exists")
                sys.exit(1)
            else:
                self.output_file = open(self.output, 'w')

    def run(self):
        if self.output != 'stdout':
            print(self.display_options())
        
        self.process_input()
        self.dump_sorted_reads()
        if self.output != 'stdout':
            print(self.display_summary())
            self.output_file.close()
    
    def process_input(self):
        """Yield a BED line from each line of a ELR input."""
        elr_in = open(self.input, 'r')
        for elr_line in elr_in:
            if elr_line[0] == '#':
                self.output_line(elr_line.rstrip())
                continue
            
            self.add_read_tuple(elr_line)
            self.linecount += 1
    
    def dump_sorted_reads(self):
        self.read_tuples.sort(reverse=True)
        while self.read_tuples:
            read_tuple = self.read_tuples.pop()
            line = '\t'.join(str(i) for i in read_tuple[:6])+'\t'+str(-read_tuple[6])
            self.output_line(line)
    
    def add_read_tuple(self, elr_line):
        l = elr_line.rstrip().split('\t')
        # parse the line as a tuple of sortable values
        read_tuple = (int(l[0]), int(l[1]), int(l[2]), l[3], l[4], int(l[5]), -float(l[6]))
        self.read_tuples.append(read_tuple)
    
    def output_line(self, line):
        """Takes a list of bed lines and writes
        them to the output stream.
        """
        if self.output_file == 'stdout':
            print(line)
        else:
            self.output_file.write('{}\n'.format(line.rstrip()))
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend sort-elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        return options_string
    
    def display_summary(self):
        summary_string = ''
        summary_string += 'Processed {} lines.\n'.format(self.linecount)
        return summary_string


if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import elr_sort_parser as parser
    args = vars(parser.parse_args())
    obj = ELRsorter(args)
    sys.exit(obj.run())
