#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy
import bookend.core.cython_utils._rnaseq_utils as ru

class GTFconverter:
    def __init__(self, args):
        """Converts each transcript of GTF or GFF3 formatted input to BED12 or ELR"""
        self.input = args['INPUT']
        self.output = args['OUT']
        self.score = args['SCORE']
        self.force = args['FORCE']
        self.gtf_parent = args['GTF_PARENT']
        self.gtf_child = args['GTF_CHILD']
        self.gff_parent = args['GFF_PARENT']
        self.gff_child = args['GFF_CHILD']
        self.parent_key_gene = args['PARENT_ATTR_GENE']
        self.child_key_gene = args['CHILD_ATTR_GENE']
        self.parent_key_transcript = args['PARENT_ATTR_TRANSCRIPT']
        self.child_key_transcript = args['CHILD_ATTR_TRANSCRIPT']
        self.color_key = args['COLOR_KEY']
        self.name_attr = args['NAME_ATTR']
        if self.color_key is None:
            self.color_code = ru.gtf_colorcode
            self.color_key = 'anno_type'
        elif args['COLOR_CODE'] is not None:
            self.color_code = self.parse_color_code(args['COLOR_CODE'])
        else:
            self.color_code = ru.gtf_colorcode
        
        self.linecount = 0
        self.outlinecount = 0

        if self.output == 'stdout':
            self.output_type = 'bed'
            self.output_file = 'stdout'
        else:
            self.output_type = self.file_extension(self.output)
            if self.output_type is None or self.output_type not in ['bed', 'bed12', 'elr']:
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
    
    @staticmethod
    def parse_color_code(filename):
        f = open(filename)
        lines = f.read()
        f.close()
        pairs = [l.split('\t') for l in lines.rstrip().split('\n')]
        return {k:v for k,v in pairs}
    
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
        for mapping_object in locus:
            self.transcript_counter += 1
            if self.output_type == 'elr':
                out_string = mapping_object.write_as_elr()
            else:
                score_column = mapping_object.attributes.get(self.score, '.')
                color = self.color_code.get(mapping_object.attributes.get(self.color_key, '.'), '80,80,80')
                out_fields = mapping_object.write_as_bed(self.dataset.chrom_array, ['.','.'], as_string=False, score_column=score_column, name_attr=self.name_attr, color=color)
                out_string = '\t'.join([str(f) for f in out_fields[:12]])
            
            self.output_line(out_string)
    
        
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
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        options_string += "  *** Parameters ***\n"
        options_string += "  Score column attribute: {}\n".format(self.score)
        input_type = self.file_extension(self.input)
        if input_type == 'gtf':
            options_string += "  Parent elements (GTF):  {}\n".format(self.gtf_parent)
            options_string += "  Child elements (GTF):   {}\n".format(self.gtf_child)
        elif input_type in ['gff','gff3']:
            options_string += "  Parent elements (GFF3): {}\n".format(self.gff_parent)
            options_string += "  Child elements (GFF3):  {}\n".format(self.gff_child)
        
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "{} loci processed ({} total transcripts).\n".format(self.locus_counter, self.transcript_counter)
        return summary

if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import gtf_to_bed_parser as parser
    args = vars(parser.parse_args())
    obj = GTFconverter(args)
    sys.exit(obj.run())
