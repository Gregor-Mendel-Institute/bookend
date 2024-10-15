#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from bookend.core.elr_combine import ELRcombiner

class ELRsorter:
    def __init__(self, args):
        """Sorts all ELR reads in a file by ascending genomic position"""
        self.input = args['INPUT']
        self.output = args['OUT']
        self.force = args['FORCE']
        self.strand_sort_values = {'+':-1, '.':0, '-':1}
        self.strand_reverse_values = {-1:'+', 0:'.', 1:'-'}
        self.read_tuples = []
        self.sortsize = 5000000
        self.linecount = 0
        self.outlinecount = 0
        self.tmpcount = 0
        self.header = ''
        if self.output == 'stdout':
            self.output_file = 'stdout'
        else:
            if self.force or not os.path.exists(self.output):
                self.output_file = open(self.output,'w')
            else:
                print("ERROR: output file already exists. Use -f/--force to overwrite.")
                sys.exit(1)

    def run(self):
        if self.output != 'stdout' and __name__ == '__main__':
            print(self.display_options())
        
        self.process_input()
        if self.tmpcount > 0:
            if self.output != 'stdout':
                print('Combining sorted reads from {} files.'.format(self.tmpcount))
            
            combine_args = {
                'INPUT':['{}.tmp{}.elr'.format(self.input,i) for i in range(self.tmpcount)],
                'OUTPUT':'stdout',
                'TEMPDIR':'{}_combinetmp'.format(self.input)
            }
            combiner = ELRcombiner(combine_args)
            for c in combiner.combine_files(combiner.input, self.output_file):pass
            
            for i in range(self.tmpcount):
                os.remove('{}.tmp{}.elr'.format(self.input,i))
        else:
            self.dump_sorted_reads()
        
        if self.output != 'stdout':
            print(self.display_summary())
            self.output_file.close()
    
    def process_input(self):
        """Yield a BED line from each line of a ELR input."""
        elr_in = open(self.input, 'r')
        for elr_line in elr_in:
            if elr_line[0] == '#':
                self.header += elr_line
                self.output_line(elr_line.rstrip())
                continue
            
            self.add_read_tuple(elr_line)
            self.linecount += 1
            if self.linecount >= self.sortsize:
                self.linecount = 0
                tmpfile = open('{}.tmp{}.elr'.format(self.input,self.tmpcount),'w')
                tmpfile.write(self.header)
                self.dump_sorted_reads(tmpfile)
                tmpfile.close()
                self.tmpcount += 1
            
        if self.tmpcount > 0:
            tmpfile = open('{}.tmp{}.elr'.format(self.input,self.tmpcount),'w')
            tmpfile.write(self.header)
            self.dump_sorted_reads(tmpfile)
            tmpfile.close()
            self.tmpcount += 1
    
    def dump_sorted_reads(self, tmpfile=None):
        """Writes sorted reads to output, collapsing
        any identical reads into a single line with increased weight"""
        self.read_tuples.sort(reverse=True)
        if len(self.read_tuples) > 0:
            read_tuple = self.read_tuples.pop()
            weights = [float(w) for w in read_tuple[6].split('|')]
            if len(weights) == 1:
                weights += [0, 0]
            
            while self.read_tuples:
                next_tuple = self.read_tuples.pop()
                if next_tuple[:6] == read_tuple[:6]:
                    new_weights = [float(w) for w in read_tuple[6].split('|')]
                    if len(new_weights) == 1:
                        new_weights += [0, 0]
                    
                    weights = [weights[0]+new_weights[0], weights[1]+new_weights[1], weights[2]+new_weights[2]]
                    continue
                else:
                    if weights[1]==0 and weights[2]==0:
                        outweight = str(round(weights[0],2))
                    else:
                        outweight = '{}|{}|{}'.format(round(weights[0],2), round(weights[1],2), round(weights[2],2))
                    
                    line = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        read_tuple[0],read_tuple[1],read_tuple[2],
                        self.strand_reverse_values[read_tuple[3]],
                        read_tuple[4],read_tuple[5],outweight
                    )
                    self.output_line(line, tmpfile)
                    self.outlinecount += 1
                    read_tuple = next_tuple
                    weights = [float(w) for w in read_tuple[6].split('|')]
                    if len(weights) == 1:
                        weights += [0, 0]
            
            if weights[1]==0 and weights[2]==0:
                outweight = str(round(weights[0],2))
            else:
                outweight = '{}|{}|{}'.format(round(weights[0],2), round(weights[1],2), round(weights[2],2))

            line = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                read_tuple[0],read_tuple[1],read_tuple[2],
                self.strand_reverse_values[read_tuple[3]],
                read_tuple[4],read_tuple[5],outweight
            )
            self.output_line(line, tmpfile)
            self.outlinecount += 1

    def add_read_tuple(self, elr_line):
        l = elr_line.rstrip().split('\t')
        # parse the line as a tuple of sortable values
        read_tuple = (int(l[0]), int(l[1]), int(l[2]), self.strand_sort_values[l[3]], l[4], int(l[5]), l[6])
        self.read_tuples.append(read_tuple)
    
    def output_line(self, line, tmpfile=None):
        """Takes a list of bed lines and writes
        them to the output stream.
        """
        if tmpfile is not None:
            tmpfile.write('{}\n'.format(line.rstrip()))
        else:
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
        summary_string += 'Processed {} lines.\n'.format((self.tmpcount*1000000)+self.linecount)
        summary_string += 'Wrote {} sorted unique reads.\n'.format(self.outlinecount)
        return summary_string


if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import elr_sort_parser as parser
    args = vars(parser.parse_args())
    obj = ELRsorter(args)
    sys.exit(obj.run())
