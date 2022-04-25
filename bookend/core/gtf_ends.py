#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy
import bookend.core.cython_utils._rnaseq_utils as ru

class GTFendwriter:
    def __init__(self, args):
        """Writes the set of nonoverlapping ends in the input annotation as BED features."""
        self.input = args['INPUT']
        self.output = args['OUT']
        self.extend = args['EXTEND']
        self.force = args['FORCE']
        self.name_attr = args['NAME_ATTR']
        self.score = args['SCORE']
        self.gtf_parent = args['GTF_PARENT']
        self.gtf_child = args['GTF_CHILD']
        self.gff_parent = args['GFF_PARENT']
        self.gff_child = args['GFF_CHILD']
        self.parent_key_gene = args['PARENT_ATTR_GENE']
        self.child_key_gene = args['CHILD_ATTR_GENE']
        self.parent_key_transcript = args['PARENT_ATTR_TRANSCRIPT']
        self.child_key_transcript = args['CHILD_ATTR_TRANSCRIPT']
        self.type = args['TYPE'].upper()
        self.source = self.input        
        self.linecount = 0
        self.outlinecount = 0
        if self.output == 'stdout':
            self.output_type = 'bed'
            self.output_file = 'stdout'
        else:
            self.output_type = self.file_extension(self.output)
            if self.output_type is None or self.output_type not in ['bed', 'bed12']:
                self.output += '.bed'
                self.output_type = 'bed'
            
            if self.force or not os.path.exists(self.output):
                self.output_file = open(self.output,'w')
            else:
                print("ERROR: output file already exists")
                sys.exit(1)
            
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        self.gtf_parent = gtf_defaults['parent_types']
        self.gtf_child = gtf_defaults['child_types']
        self.gff_parent = gff_defaults['parent_types']
        self.gff_child = gff_defaults['child_types']
        self.dataset = ru.AnnotationDataset(
            annotation_files=[],
            reference=self.input, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        self.dataset.source_array = [self.source]
        self.generator = self.dataset.generator
        self.locus_counter = 0
        self.transcript_counter = 0
    
    def run(self):
        if self.output != 'stdout':
            print(self.display_options())
                
        for locus in self.generator:
            self.locus_counter += 1
            self.process_locus(locus)
        
        if self.output != 'stdout':
            print(self.display_summary())
        
    @staticmethod
    def file_extension(filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        else:
            extension = split_name[-1].lower()
            return extension
    
    def make_config_dicts(self):
        """Converts commandline input into three config dicts
        to pass to the AnnotationDataset."""
        config_defaults = copy.copy(ru.config_defaults)
        gtf_defaults = copy.copy(ru.gtf_defaults)
        gff_defaults = copy.copy(ru.gff_defaults)
        config_defaults['min_reps'] = 1
        config_defaults['cap_percent'] = 0
        config_defaults['verbose'] = False
        if self.gtf_parent: gtf_defaults['parent_types'] = set(self.gtf_parent)
        if self.gtf_child: gtf_defaults['child_types'] = set(self.gtf_child)
        if self.gff_parent: gff_defaults['parent_types'] = set(self.gff_parent)
        if self.gff_child: gff_defaults['child_types'] = set(self.gff_child)
        if self.parent_key_transcript is not None:
            gtf_defaults['parent_key_transcript'] += self.parent_key_transcript
            gff_defaults['parent_key_transcript'] += self.parent_key_transcript
        
        if self.child_key_transcript is not None:
            gtf_defaults['child_key_transcript'] += self.child_key_transcript
            gff_defaults['child_key_transcript'] += self.child_key_transcript
        
        if self.parent_key_gene is not None:
            gtf_defaults['parent_key_gene'] = self.parent_key_gene
            gff_defaults['parent_key_gene'] = self.parent_key_gene
        
        if self.child_key_gene is not None:
            gtf_defaults['child_key_gene'] = self.child_key_gene
            gff_defaults['child_key_gene'] = self.child_key_gene    
        
        return config_defaults, gtf_defaults, gff_defaults
    
    def process_locus(self, locus):
        """Given a chunk of transcripts from an AnnotationDataset, print all as BED12"""
        # start/end are tuples (leftmost, rightmost, peak, name)
        SP, SM, EP, EM = [], [], [], []
        for transcript in locus:
            chrom = self.dataset.chrom_array[transcript.chrom]
            self.transcript_counter += 1
            name = transcript.attributes.get(self.name_attr, None)
            if name is None:
                name = self.transcript_counter
            
            if transcript.strand == 1:
                startname = str(name)
                endname = str(name)
                if self.extend is None:
                    pos = transcript.left()
                    start = (int(transcript.attributes.get('S.left',pos)), int(transcript.attributes.get('S.right',pos)), pos, startname) 
                    pos = transcript.right()
                    end = (int(transcript.attributes.get('E.left',pos)), int(transcript.attributes.get('E.right',pos)), pos, endname)
                else:
                    pos = transcript.left()
                    start = (pos-self.extend, pos+self.extend, pos, startname)
                    pos = transcript.right()
                    end = (pos-self.extend, pos+self.extend, pos, endname)
                
                SP.append(start)
                EP.append(end)
            elif transcript.strand == -1:
                startname = str(name)
                endname = str(name)
                if self.extend is None:
                    pos = transcript.right()
                    start = (int(transcript.attributes.get('S.left',pos)), int(transcript.attributes.get('S.right',pos)), pos, startname) 
                    pos = transcript.left()
                    end = (int(transcript.attributes.get('E.left',pos)), int(transcript.attributes.get('E.right',pos)), pos, endname)
                else:
                    pos = transcript.right()
                    start = (pos-self.extend, pos+self.extend, pos, startname)
                    pos = transcript.left()
                    end = (pos-self.extend, pos+self.extend, pos, endname)
                
                SM.append(start)
                EM.append(end)
            
        SPout, SMout, EPout, EMout = [], [], [], []
        if self.type in ['ALL','S']:
            SPout = self.merge_ends(sorted(SP), 'SP')
            SMout = self.merge_ends(sorted(SM), 'SM')

        if self.type in ['ALL','E']:
            EPout = self.merge_ends(sorted(EP), 'EP')
            EMout = self.merge_ends(sorted(EM), 'EM')
        
        for bed in sorted([(l,r,'SP:'+name,peak,'+') for l,r,peak,name in SPout] +
            [(l,r,'SM:'+name,peak,'-') for l,r,peak,name in SMout] +
            [(l,r,'EP:'+name,peak,'+') for l,r,peak,name in EPout] +
            [(l,r,'EM:'+name,peak,'-') for l,r,peak,name in EMout]):
            self.linecount += 1
            self.output_line('{}\t{}\t{}\t{}\t{}\t{}'.format(
                chrom, bed[0], bed[1], bed[2], bed[3], bed[4]
            ))
    
    def merge_ends(self, end_tuples, endtype):
        """Given a list of (left, right, peak, name) tuples,
        merge into a unique set of nonoverlapping features."""
        peakfunction = {'SP':min, 'SM':max, 'EP':max, 'EM':min}[endtype]
        mergelist = []
        if len(end_tuples) > 0:
            last_tuple = end_tuples[0]
            for i in range(1,len(end_tuples)):
                next_tuple = end_tuples[i]
                if next_tuple[0] > last_tuple[1]:
                    mergelist.append(last_tuple)
                    last_tuple = next_tuple
                else: # merge if overlap
                    last_tuple = (
                        min(last_tuple[0],next_tuple[0]),
                        max(last_tuple[1], next_tuple[1]),
                        peakfunction(last_tuple[2],next_tuple[2]),
                        ','.join([last_tuple[3], next_tuple[3]])
                    )
            
            mergelist.append(last_tuple)
        
        return mergelist
    
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
        options_string = "\n/| bookend gtf-to-bed |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:     {}\n".format(self.input)
        options_string += "  Output file:    {}\n".format(self.output)
        options_string += "  *** Parameters ***\n"
        options_string += "  Name attribute: {}\n".format(self.name_attr)
        options_string += "  End type:       {}\n".format(self.type)
        input_type = self.file_extension(self.input)
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "{} transcripts processed ({} unique ends).\n".format(self.transcript_counter, self.linecount)
        return summary

if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import gtf_ends_parser as parser
    args = vars(parser.parse_args())
    obj = GTFendwriter(args)
    sys.exit(obj.run())
