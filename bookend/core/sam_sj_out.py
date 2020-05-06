#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pysam
if __name__ == '__main__':
    sys.path.append('../../bookend')

import bookend.core.cython_utils._fasta_utils as fu
import bookend.core.cython_utils._rnaseq_utils as ru
from bookend.core.sj_merge import SJobject

class SAMtoSJconverter:
    def __init__(self, args):
        """Generates an SJ.out.tab file with splice junction information from a SAM file."""
        self.junction_types = {'GTAG':1, 'CTAC':2, 'GCAG':3, 'CTGC':4, 'ATAC':5, 'GTAT':6}    
        self.fasta = args['FASTA']
        self.format = args['FORMAT']
        self.filter = args['FILTER']
        self.input = args['INPUT']
        self.sj_dict = {}
        self.sam_in = pysam.AlignmentFile(self.input)
        self.dataset = ru.RNAseqDataset(
            chrom_array=self.sam_in.header.references, 
            chrom_lengths=list(self.sam_in.header.lengths),
            source_array=['SAM'],
            genome_fasta=self.fasta
        )
        
        self.generator = self.generate_bam_entries()
    
    def run(self):
        print(self.display_options())
        for entry in self.generator:
            self.add_entry_to_sj_dict(entry)
        
        self.write_sj_dict_to_file()
        print(self.display_summary())
    
    def write_sj_dict_to_file(self):
        outfile_name = self.input+'.SJ.out.tab'
        output_file = open(outfile_name, 'w')
        for sj in sorted(list(self.sj_dict.items())):
            output_file.write(sj.as_string(self.format))
        
        output_file.close()
    
    def generate_bam_entries(self):
        """Group all BAM lines with the same read ID
        and yield them as a list of pysam objects"""
        current_ID = None
        bam_lines = []
        for line in self.sam_in:
            ID = line.query_name
            if ID == current_ID or current_ID is None:
                bam_lines.append(line)
                current_ID = ID
            else:
                yield bam_lines
                bam_lines = [line]
                current_ID = ID
        
        yield bam_lines
    
    def get_junction_type(self, chrom, left, right):
        """ Returns the sequence motif ID for a splice junction based
        on the flanking genomic sequence."""
        flanking_sequence = ru.get_flank(self.dataset.genome, chrom, left-1, '+', 'E', 2) + ru.get_flank(self.dataset.genome, chrom, right, '+', 'S', 2)
        flanking_sequence = flanking_sequence.upper()
        jtype = self.junction_types.get(flanking_sequence,0)
        return jtype
    
    def add_entry_to_sj_dict(self, entry):
        """Extract all splice junctions from the read(s) of a SAM entry"""
        try:
            Nmap = int(entry[0].get_tag('NH'))
        except KeyError:
            if entry[0].is_paired:
                Nmap = int(len(entry)*0.5)
            else:
                Nmap = len(entry)
        
        if Nmap > 1: # multimapper
            unique, multi = 0, 1
        elif Nmap == 1: # unique mapper
            unique, multi = 1, 0
        
        self.dataset.add_read_from_BAM(entry)
        read = self.dataset.read_list.pop()
        chrom = self.dataset.chrom_array[read.chrom]
        strand = 2 if read.strand == -1 else read.strand
        junctions = read.junctions()
        for i in range(len(junctions)):
            l,r = junctions[i]
            minoverhang = min([l-read.ranges[i][0], read.ranges[i+1][1]-r])
            motif = self.get_junction_type(chrom, l, r)
            if self.filter and motif not in [1,2]:
                continue
            
            SJ_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                chrom, l+1, r, strand, motif, 0, unique, multi, minoverhang
            )
            sj = SJobject(SJ_line, 'star')
            if sj.hash in self.sj_dict:
                self.sj_dict[sj.hash].merge(sj)
            else:
                self.sj_dict[sj.hash] = sj
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend index-fasta |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:         {}\n".format(self.input)
        options_string += "  Genome FASTA:       {}\n".format(self.fasta)
        options_string += "  *** Parameters ***\n"
        options_string += "  Output format:      {}\n".format(self.format)
        options_string += "  Canonical SJ only:  {}\n".format(self.bridge_min)
        return options_string        
    
    def display_summary(self):
        summary_string = ''
        summary_string += 'Splice junctions written to {}.\n'.format(self.input+'.SJ.out.tab')
        
        return summary_string

if __name__ == '__main__':
    from argument_parsers import sam_sj_parser as parser
    args = vars(parser.parse_args())
    obj = SAMtoSJconverter(args)
    sys.exit(obj.run())