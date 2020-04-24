#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import bookend.core.cython_utils._fasta_utils as fu
if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import fasta_index_parser as parser

class Indexer:
    def __init__(self, args):
        """Generates a set of three index files for an input genome FASTA file."""
        self.fasta = args['FASTA']
        self.split_on = args['SPLIT_ON']
        self.bridge_min = args['MINLEN']
        self.bridge_max = args['MAXLEN']
        self.index = self.fasta+'.fai'
        self.lengths = self.fasta+'.lengths.tsv'
        self.softbridges = self.fasta+'.softbridges.bed'
        self.index_exists = os.path.exists(self.index)
        self.lengths_exists = os.path.exists(self.lengths)
        self.softbridges_exists = os.path.exists(self.softbridges)
    
    def run(self):
        print(self.display_options())
        if not (self.index_exists and self.lengths_exists and self.softbridges_exists):
            self.genome, self.index_string = fu.import_genome(self.fasta, self.split_on, True, self.index_exists)
            if not self.index_exists:
                self.write_index()
            
            if not self.lengths_exists:
                self.write_lengths()
            
            if not self.softbridges_exists:
                self.write_softbridges()
            
        print(self.display_summary())
        return 0
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend index-fasta |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:          {}\n".format(self.fasta)
        options_string += "  Output index:        {}\n".format(self.index)
        options_string += "  Output length table: {}\n".format(self.lengths)
        options_string += "  Output softbridges:  {}\n".format(self.softbridges)
        options_string += "  *** Parameters ***\n"
        options_string += "  Name split character (--split):       {}\n".format(self.split_on)
        options_string += "  Min softbridge length (--bridge_min): {}\n".format(self.bridge_min)
        options_string += "  Max softbridge length (--bridge_max): {}\n".format(self.bridge_max)
        return options_string
    
    def display_summary(self):
        summary_string = ''
        if self.index_exists:
            summary_string += 'Index file ({}) already exists.\n'.format(self.index)
        else:
            summary_string += 'Index written to ({}).\n'.format(self.index)
        
        if self.lengths_exists:
            summary_string += 'Length table ({}) already exists.\n'.format(self.lengths)
        else:
            summary_string += 'Length table written to ({}).\n'.format(self.lengths)
        
        if self.softbridges_exists:
            summary_string += 'Softbridge file ({}) already exists.\n'.format(self.softbridges)
        else:
            summary_string += 'Softbridges written to ({}).\n'.format(self.softbridges)
        
        return summary_string
    
    def write_index(self):
        index_file = open(self.index, 'w')
        index_file.write(self.index_string)
        index_file.close()
    
    def write_lengths(self):
        lengths_string = ''
        chromosomes = {}
        for chromname, chromstring in self.genome.items():
            chromosomes[chromname] = len(chromstring)
        
        for chromname, chromlength in sorted(list(chromosomes.items())):
            lengths_string += '{}\t{}\n'.format(chromname, chromlength)
        
        self.lengths_file = open(self.lengths, 'w')
        self.lengths_file.write(lengths_string)
        self.lengths_file.close()
    
    def write_softbridges(self):
        self.softbridges_file = open(self.softbridges, 'w')
        sb_gen = fu.generate_softbridges(self.genome, self.bridge_min, self.bridge_max)
        for bridge in sb_gen:
            self.softbridges_file.write(bridge)
        
        self.softbridges_file.close()

if __name__ == '__main__':
    args = vars(parser.parse_args())
    obj = Indexer(args)
    sys.exit(obj.run())