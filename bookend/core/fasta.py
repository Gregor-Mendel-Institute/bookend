#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import copy
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._fasta_utils as fu
if __name__ == '__main__':
    sys.path.append('../../bookend')

class FastaWriter:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.args = args
        self.output = self.args['OUT']
        self.genome = self.args['GENOME']
        self.input = self.args['INPUT']
        self.orf = self.args['ORF']
        self.allow_unstranded = self.args['UNSTRANDED']
        if len(self.input) == 0:
            parser.print_help()
            sys.exit(0)
        
        if self.file_extension(self.output) not in ['fa','fasta']: # Check for valid file extension on output name
            self.output_type = 'fasta'
            self.output = self.output + '.fasta'
        else:
            self.output_file = open(self.output,'w')
            self.output_type = self.file_extension(self.output)
            self.output_file = open(self.output,'w')
        
        print(self.display_options())
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        self.dataset = ru.AnnotationDataset(
            annotation_files=self.input, 
            reference=None, 
            genome_fasta=self.genome, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        self.dataset.source_array = ['', 'bookend']
        self.locus_counter = 0
        self.new_gene_counter = 0
        self.input_transcripts = 0
        self.transcript_counter = 0
        self.updated_transcript_counter = 0
    
    def make_config_dicts(self):
        """Converts commandline input into three config dicts
        to pass to the AnnotationDataset."""
        config_defaults = copy.copy(ru.config_defaults)
        gtf_defaults = copy.copy(ru.gtf_defaults)
        gff_defaults = copy.copy(ru.gff_defaults)
        return config_defaults, gtf_defaults, gff_defaults
    
    def output_transcripts(self, transcript):
        """Writes the transcript FASTA entry corresponding to an RNAseqMapping object."""
        sequence = self.dataset.get_transcript_fasta(transcript)
        if self.orf:
            orf = fu.longest_orf(sequence)[0]
            if transcript.strand == 0:
                rev_orf = fu.longest_orf(fu.rc(sequence))[0]
                if len(rev_orf) > len(orf):
                    orf = rev_orf
            
            sequence = orf
        
        if len(sequence) == 0:
            sequence = '-'
        
        self.output_file.write('>{}\n{}\n'.format(
            transcript.attributes['transcript_id'],
            sequence
        ))
    
    def process_entry(self, chunk):
        self.locus_counter += 1
        for transcript in chunk:
            if not self.allow_unstranded and transcript.strand == 0:
                continue
            
            self.transcript_counter += 1
            self.output_transcripts(transcript)
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend fasta |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:\n\t{}\n".format('\n\t'.join(self.input))
        options_string += "  Reference genome (--genome):\n\t{}\n".format(self.genome)
        options_string += "  Output file (-o):\n\t{}\n".format(self.output)
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "{} loci processed ({} total input transcripts).\n".format(self.locus_counter, self.input_transcripts)
        summary += "{} transcripts written.\n".format(self.transcript_counter)
        return summary
    
    def file_extension(self, filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        else:
            extension = split_name[-1].lower()
            return extension
    
    def input_is_valid(self, filename):
        """Boolean if the file is a format that Assembler can parse."""
        if self.file_extension(filename) in ['bed','bed12','elr','gtf','gff3','gff']:
            return True
        else:
            return False
    
    def navigate_to(self, gene_id):
        generator = self.dataset.generate_loci()
        for chunk in generator:
            if gene_id in [r.attributes['gene_id'] for r in chunk]:
                return copy.deepcopy(chunk)

    def run(self):
        """Executes end labeling on all reads."""
        for chunk in self.dataset.generator:
            if len(chunk) > 0:
                self.process_entry(chunk)
        
        self.output_file.close()
        print(self.display_summary())


if __name__ == '__main__':
    from argument_parsers import merge_parser as parser
    args = vars(parser.parse_args())
    obj = FastaWriter(args)
    sys.exit(obj.run())

