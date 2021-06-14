#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import copy
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._assembly_utils as au
import bookend.core.cython_utils._fasta_utils as fu
import numpy as np
from collections import Counter, namedtuple
from math import ceil
if __name__ == '__main__':
    sys.path.append('../../bookend')


class AssemblyClassifier:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.match_data = namedtuple('match_data', 'matchtype transcript gene exonoverlap reflen tlen ref diff5p diff3p')
        self.args = args
        self.output = self.args['OUT']
        self.end_buffer = self.args['END_BUFFER']
        self.input = self.args['INPUT']
        self.verbose = self.args['VERBOSE']
        self.reference = self.args['REFERENCE']
        if len(self.input) == 0 and self.reference is None:
            parser.print_help()
            sys.exit(0)
        
        self.gene_attr = self.args['PARENT_ATTR_GENE']
        self.gene_attr_child = self.args['CHILD_ATTR_GENE']
        self.gtf_parent = self.args['GFF_PARENT']
        self.gtf_child = self.args['GFF_CHILD']
        self.gff_parent = self.args['GFF_PARENT']
        self.gff_child = self.args['GFF_CHILD']
        self.refid_parent = self.args['PARENT_ATTR_TRANSCRIPT']
        self.refid_child = self.args['CHILD_ATTR_TRANSCRIPT']
        self.gene_delim = self.args['GENE_DELIM']
        if self.input_is_valid(self.output): # Check for valid file extension on output name
            self.output_type = self.file_extension(self.output)
            self.output_file = open(self.output,'w')
        else:
            self.output_type = 'tsv'
            self.output = self.output + '.tsv'
            self.output_file = open(self.output,'w')
        
        print(self.display_options())
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        self.dataset = ru.AnnotationDataset(
            annotation_files=self.input, 
            reference=self.reference, 
            genome_fasta=None, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        self.output_file.write('assembly_id\tclassification\tref_match\tref_gene\tassembly_len\tref_len\toverlap_len\tdiff5p\tdiff3p\tcov\tS.reads\tS.capped\tE.reads\n')
        self.dataset.source_array = ['reference', 'assembly']
        self.generator = self.dataset.generator
        self.locus_counter = 0
        self.new_gene_counter = 0
        self.ref_transcript_count = 0
        self.input_transcript_count = 0
        self.transcript_counter = 0
        self.updated_transcript_counter = 0
        self.match_types = [
            'intergenic', # 0 (lowest classification) no ref match
            'ambiguous',  # 1 overlapping with no strand information
            'antisense',  # 2 only overlaps a ref in antisense
            'intronic',   # 3 fully contained in a ref intron (sense)
            'isoform',    # 4 overlaps, incompatible exon chain
            'fragment',   # 5 compatible with, but fewer exons than, a ref
            'fusion',     # 6 shares exons with 2 or more ref genes
            'exon_match', # 7 shares entire exon chain, but not ends
            'full_match'  # 8 shares entire exon chan and ends
        ]
    
    def make_config_dicts(self):
        """Converts commandline input into three config dicts
        to pass to the AnnotationDataset."""
        config_defaults = copy.copy(ru.config_defaults)
        gtf_defaults = copy.copy(ru.gtf_defaults)
        gff_defaults = copy.copy(ru.gff_defaults)
        config_defaults['min_reps'] = 0
        config_defaults['cap_percent'] = 0
        config_defaults['verbose'] = self.verbose
        config_defaults['gene_delim'] = self.gene_delim
        if self.gtf_parent: gtf_defaults['parent_types'] = set(self.gtf_parent)
        if self.gtf_child: gtf_defaults['child_types'] = set(self.gtf_child)
        if self.gff_parent: gff_defaults['parent_types'] = set(self.gff_parent)
        if self.gff_child: gff_defaults['child_types'] = set(self.gff_child)
        if self.gene_attr is None:
            self.gene_attr = gtf_defaults['parent_key_gene']
        else:
            gff_defaults['parent_key_gene'] = self.gene_attr
            gtf_defaults['parent_key_gene'] = self.gene_attr
        
        if self.gene_attr_child is not None:
            gff_defaults['child_key_gene'] = self.gene_attr_child
        
        if self.refid_parent is not None:
            gff_defaults['parent_key_transcript'] += self.refid_parent
            gtf_defaults['parent_key_transcript'] += self.refid_parent
        
        if self.refid_child is not None:
            gff_defaults['child_key_transcript'] += self.refid_child
            gtf_defaults['child_key_transcript'] += self.refid_child
        
        return config_defaults, gtf_defaults, gff_defaults
    
    def compare_transcripts(self, input_transcripts, reference_transcripts):
        """Iterate over the set of input transcripts (highest TPM first)
        and classify each input relative to the reference according to the
        definition in calculate_match_type()
        """
        for transcript in input_transcripts:
            match_data = self.calculate_match_type(transcript, reference_transcripts)
            classification = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                transcript.attributes['transcript_id'],
                self.match_types[match_data.matchtype],
                match_data.transcript,
                match_data.gene,
                match_data.tlen,
                match_data.reflen,
                match_data.exonoverlap,
                match_data.diff5p,
                match_data.diff3p,
                round(float(transcript.attributes.get('cov', 0)),1),
                round(float(transcript.attributes.get('S.reads', 0)),1),
                round(float(transcript.attributes.get('S.capped', 0)),1),
                round(float(transcript.attributes.get('E.reads', 0)),1)
            )
            self.output_file.write(classification)
            if self.verbose:
                print(classification.rstrip())
        
        return
    
    def update_match(self, old_match, new_match):
        """Given an existing ref match, decide how the new match
        updates the old match. Can replace the old match, be skipped, 
        or turn the match into a fusion the transcript matches two genes.
            'matchtype transcript gene exonoverlap reflen tlen'
        """
        if old_match.matchtype == 8:
            if new_match.matchtype == 8 and new_match.exonoverlap > old_match.exonoverlap:
                return new_match
            else:
                return old_match
        elif new_match.matchtype == 8:
            return new_match
        elif new_match.gene != old_match.gene and old_match.gene != 'NA': # Evaluate if a fusion
            if new_match.exonoverlap >= .9*new_match.reflen or new_match.matchtype == 4:
                if old_match.exonoverlap >= .9*old_match.reflen or old_match.matchtype == 4:
                    if old_match.ref.shared_bases(new_match.ref) < .9*min([old_match.reflen, new_match.reflen]):
                        fused_match = self.match_data(
                            6, 
                            '{},{}'.format(old_match.transcript, new_match.transcript), 
                            '{},{}'.format(old_match.gene, new_match.gene),
                            max(old_match.exonoverlap, new_match.exonoverlap),
                            max(old_match.reflen, new_match.reflen),
                            old_match.tlen,
                            old_match.ref if old_match.exonoverlap > new_match.exonoverlap else new_match.ref,
                            old_match.diff5p if abs(old_match.diff5p) <= abs(new_match.diff5p) else new_match.diff5p,
                            old_match.diff3p if abs(old_match.diff3p) <= abs(new_match.diff3p) else new_match.diff3p
                        )
                        return fused_match
        
        if new_match.matchtype > old_match.matchtype:
            return new_match
        elif new_match.matchtype == old_match.matchtype and new_match.exonoverlap > old_match.exonoverlap:
            return new_match
        else:
            return old_match
    
    def calculate_match_type(self, transcript, reference_transcripts):
        """Finds the best match in the list of reference transcripts
        and classifies it according to one of the self.match_types:
            'intergenic', # 0 (lowest classification) no ref match
            'ambiguous',  # 1 overlapping with no strand information
            'antisense',  # 2 only overlaps a ref in antisense
            'intronic',   # 3 fully contained in a ref intron (sense)
            'isoform',    # 4 overlaps, incompatible exon chain
            'truncation', # 5 compatible with, but fewer exons than, a ref
            'fusion',     # 6 shares exons with 2 or more ref genes
            'exon_match', # 7 shares entire exon chain, but not ends
            'full_match'  # 8 shares entire exon chan and ends
        """
        tlen = transcript.get_length()
        best_match = self.match_data(0, 'NA', 'NA', 0, 0, tlen, None, 'NA', 'NA')
        for ref in reference_transcripts:
            match_type = 0
            if not transcript.overlaps(ref):
                continue
            
            reflen = ref.get_length()
            shared_bases = transcript.shared_bases(ref)
            left_diff = transcript.ranges[0][0]-ref.ranges[0][0]
            right_diff = transcript.ranges[-1][1]-ref.ranges[-1][1]
            diff5p = left_diff if transcript.strand == 1 else right_diff
            diff3p = right_diff if transcript.strand == 1 else left_diff
            if transcript.splice_match(ref, ignore_ends=True):
                if abs(left_diff) <= self.end_buffer and abs(right_diff) <= self.end_buffer:
                    match_type = 8 # full_match
                else:
                    match_type = 7 # exon_match
            elif transcript.is_compatible(ref, ignore_ends=True, ignore_source=True): # Truncation or extension
                if reflen > tlen:
                    match_type = 5 # truncation
                else:
                    match_type = 4 # isoform (possibly a fusion of two genes)
            else: # Incompatible but overlapping
                if transcript.strand == 0:
                    match_type = 1 # ambiguous
                elif transcript.strand != ref.strand:
                    match_type = 2 # antisense
                elif shared_bases == 0:
                    match_type = 3
                else: # At least some shared, sense, exonic sequence
                    match_type = 4 # isoform
            
            if transcript.strand == 0 and match_type in [8, 7, 4]:
                match_type = 1
            
            new_match = self.match_data(match_type, ref.attributes['transcript_id'], ref.attributes[self.gene_attr], shared_bases, reflen, tlen, ref, diff5p, diff3p)
            best_match = self.update_match(best_match, new_match)
        
        return best_match
    
    def process_entry(self, chunk):
        self.locus_counter += 1
        ref_transcripts = [t for t in chunk if t.is_reference]
        self.ref_transcript_count += len(ref_transcripts)
        input_transcripts = [t for t in chunk if not t.is_reference]
        self.input_transcript_count += len(input_transcripts)
        if len(input_transcripts) > 0: # Work needs to be done on non-reference transcripts
            self.compare_transcripts(input_transcripts, ref_transcripts)
        
        return
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend classify |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input files:\n\t{}\n".format('\n\t'.join(self.input))
        options_string += "  Reference file (-r):\n\t{}\n".format(self.reference)
        options_string += "  Output file (-o):\n\t{}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Cluster distance for ends (--end_buffer):   {}\n".format(self.end_buffer)
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "{} loci processed ({} reference transcripts, {} input transcripts).\n".format(self.locus_counter, self.ref_transcript_count, self.input_transcript_count)
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
        if self.file_extension(filename) in ['tsv','txt']:
            return True
        else:
            return False
    
    def navigate_to(self, gene_id):
        generator = self.dataset.generate_loci()
        for chunk in generator:
            if gene_id in [r.attributes[self.gene_attr] for r in chunk]:
                return copy.deepcopy(chunk)

    def run(self):
        """Executes end labeling on all reads."""
        for chunk in self.generator:
            if len(chunk)>0:
                self.process_entry(chunk)
        
        self.output_file.close()
        print(self.display_summary())


if __name__ == '__main__':
    from argument_parsers import classify_parser as parser
    args = vars(parser.parse_args())
    obj = AssemblyClassifier(args)
    sys.exit(obj.run())

