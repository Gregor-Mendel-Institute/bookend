#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import resource
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._fasta_utils as fu
from bookend.core.cython_utils._pq import IndexMinPQ
if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import combine_parser as parser

class ELRcombiner:
    def __init__(self):
        """Leaves together ELR files in sort order"""
        self.input = args['INPUT']
        self.output = args['OUTPUT']
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            self.output_file = open(self.output, 'w')
        
        self.linecount = 0
        self.readcount = 0
        self.dataset = ru.RNAseqDataset()

    def get_header(self, file):
        """From an open connection to an ELR file, 
        extract the header information and store it as
        a master lookup table for output chrom and source.
        Returns the first line without a header and keeps
        the file open at the existing buffer position."""
        header = {'chrom':{}, 'source':{}}
        file.seek(0)
        line = file.readline().rstrip()
        if len(line) == 0:
            return header, line
        
        while line[0] == '#':
            header_line = line.split(' ')
            num, char = header_line[-2:]
            if header_line[0] == '#S':
                header['source'][num] = char
            elif header_line[0] == '#C':
                header['chrom'][num] = char
            
            line = file.readline().rstrip()
        
        return header, line

    def read_to_sortable_tuple(self, line, index):
        """Convert an ELR line into an RNAseqMapping object
        with shared chrom and source indices."""
        split_line = line.split('\t')
        # Swap index for the merged index value
        chrom = int(dataset.chrom_dict[file_headers[index]['chrom'][split_line[0]]])
        start = int(split_line[1])
        length = int(split_line[2])
        strand = strand_sort_values[split_line[3]]
        elcigar = split_line[4]
        source = int(dataset.source_dict[file_headers[index]['source'][split_line[5]]])
        weight = -float(split_line[6])
        return (chrom, start, length, strand, elcigar, source, weight)

    def sortable_tuple_to_read(self, sortable_tuple):
        l = list(sortable_tuple)
        l[3] = strand_reverse_values[l[3]]    
        l[6] = -l[6]
        return '\t'.join([str(i) for i in l])
    
    def run(self):
        if self.output != 'stdout':
            self.display_options()
        
        self.combine_files()
        if self.output != 'stdout':
            self.display_summary()
    
    def combine_files(self):
        if args.INPUT:
            if not all([i[-3:].lower()=='elr' for i in args.INPUT]):
                print("\nERROR: all input files must be ELR format.")
                parser.print_help()
                sys.exit(1)
        
            filenames = args.INPUT
            #TODO: Check open file resource limits!
            files = [open(f) for f in filenames]
        else:
            print("\nERROR: requires ELR file as input.")
            parser.print_help()
            sys.exit(1)
        
        strand_sort_values = {'+':-1, '.':0, '-':1}
        strand_reverse_values = {-1:'+', 0:'.', 1:'-'}
        number_of_files = len(filenames)
        file_number = 0
        file_headers = ['']*number_of_files
        current_lines = ['']*number_of_files
        for i in range(number_of_files):
            file_headers[i], current_lines[i] = self.get_header(files[i])
        
        set_of_chroms = set()
        set_of_sources = set()
        for h in file_headers:
            set_of_chroms.update(h['chrom'].values())
            set_of_sources.update(h['source'].values())
        
        merged_chroms = sorted(list(set_of_chroms))
        merged_sources = sorted(list(set_of_sources))
        dataset = ru.RNAseqDataset(chrom_array=merged_chroms, source_array=merged_sources)
        # Initialize the priority queue with one ELR read object from each file
        PQ = IndexMinPQ(number_of_files)
        finished_files = 0
        for index, line in enumerate(current_lines):
            if line: # Populate the queue with the first line
                item = self.read_to_sortable_tuple(line, index)
                PQ.insert(index, item)
            else: # The file was empty
                files[index].close()
                finished_files += 1
        
        for h in dataset.dump_header():
            print(h)
        
        while finished_files < number_of_files: # Keep going until every line of every file is processed
            index, item = PQ.pop(True)
            print(self.sortable_tuple_to_read(item))
            next_line = files[index].readline().rstrip()
            if next_line:
                next_item = self.read_to_sortable_tuple(next_line, index)
                PQ.insert(index, next_item)
            else: # No more lines in this file
                files[index].close()
                finished_files += 1
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend bed-to-elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:     {}\n".format(self.input)
        options_string += "  Output file:    {}\n".format(self.output)
        options_string += "  Temp directory: {}\n".format(self.header)
        return options_string
    
    def display_summary(self):
        summary_string = ''
        summary_string += 'Wrote {} lines ({} total reads) from {} files.\n'.format(self.linecount, round(self.readcount,2), len(self.input))
        return summary_string

####################
# PARSE INPUT FILE #
####################
if __name__ == '__main__':
    args = vars(parser.parse_args())
    obj = ELRcombiner(args)
    sys.exit(obj.run())
