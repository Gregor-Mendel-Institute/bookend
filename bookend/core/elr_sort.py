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
        self.outlinecount = 0
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            if self.force or not os.path.exists(self.output):
                self.output_file = open(self.output,'w')
            else:
                print("ERROR: output file already exists")
                sys.exit(1)

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
        """Writes sorted reads to output, collapsing
        any identical reads into a single line with increased weight"""
        self.read_tuples.sort(reverse=True)
        if len(self.read_tuples) > 0:
            read_tuple = self.read_tuples.pop()
            weight = -read_tuple[6]
            while self.read_tuples:
                next_tuple = self.read_tuples.pop()
                if next_tuple[:6] == read_tuple[:6]:
                    weight += -next_tuple[6]
                    continue
                else:
                    line = '\t'.join(str(i) for i in read_tuple[:6])+'\t'+str(weight)
                    self.output_line(line)
                    self.outlinecount += 1
                    read_tuple = next_tuple
                    weight = -read_tuple[6]
            
            line = '\t'.join(str(i) for i in read_tuple[:6])+'\t'+str(weight)
            self.output_line(line)
            self.outlinecount += 1

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
        summary_string += 'Wrote {} sorted unique reads.\n'.format(self.outlinecount)
        return summary_string


if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import elr_sort_parser as parser
    args = vars(parser.parse_args())
    obj = ELRsorter(args)
    sys.exit(obj.run())
