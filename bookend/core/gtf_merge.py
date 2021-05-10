#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import copy
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._assembly_utils as au
import bookend.core.cython_utils._fasta_utils as fu
import numpy as np
from collections import Counter
from math import ceil
if __name__ == '__main__':
    sys.path.append('../../bookend')

class AnnotationMerger:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.args = args
        self.output = self.args['OUT']
        self.fasta_out = self.args['FASTA_OUT']
        self.orf_out = self.args['ORF_OUT']
        self.genome = self.args['GENOME']
        self.end_cluster = self.args['END_CLUSTER']
        self.min_reps = self.args['MIN_REPS']
        self.minlen = self.args['MINLEN']
        self.cap_percent = self.args['CAP_PERCENT']
        self.input = self.args['INPUT']
        self.confidence_threshold = self.args['CONFIDENCE']
        self.confidence =  ceil(self.confidence_threshold * len(self.input))
        self.verbose = self.args['VERBOSE']
        self.reference = self.args['REFERENCE']
        self.refname = self.args['REFNAME']
        if self.refname is None and self.reference is not None:
            self.refname = '.'.join(self.reference.split('/')[-1].split('.')[:-1])
        
        if len(self.input) == 0 and self.reference is None:
            parser.print_help()
            sys.exit(0)
        
        self.gtf_parent = self.args['GTF_PARENT']
        self.gtf_child = self.args['GTF_CHILD']
        self.gff_parent = self.args['GFF_PARENT']
        self.gff_child = self.args['GFF_CHILD']
        self.refid_parent = self.args['REF_ID_PARENT']
        self.refid_child = self.args['REF_ID_CHILD']

        if self.input_is_valid(self.output): # Check for valid file extension on output name
            self.output_type = self.file_extension(self.output)
            self.output_file = open(self.output,'w')
        else:
            self.output_type = 'gtf'
            self.output = self.output + '.gtf'
            self.output_file = open(self.output,'w')
        
        if self.fasta_out:
            self.output_file_fasta = open(self.fasta_out, 'w')
        
        if self.orf_out:
            self.output_file_orf = open(self.orf_out, 'w')
        
        print(self.display_options())
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        self.dataset = ru.AnnotationDataset(
            annotation_files=self.input, 
            reference=self.reference, 
            genome_fasta=self.genome, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        self.dataset.source_array = [self.refname, 'bookend']
        self.generator = self.dataset.generator
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
        config_defaults['min_reps'] = self.min_reps
        config_defaults['cap_percent'] = self.cap_percent
        config_defaults['verbose'] = self.verbose
        if self.gtf_parent: gtf_defaults['parent_types'] = set(self.gtf_parent)
        if self.gtf_child: gtf_defaults['child_types'] = set(self.gtf_child)
        if self.gff_parent: gff_defaults['parent_types'] = set(self.gff_parent)
        if self.gff_child: gff_defaults['child_types'] = set(self.gff_child)
        if len(self.refid_parent)>0:
            gff_defaults['parent_key_transcript'] += self.refid_parent
            gtf_defaults['parent_key_transcript'] += self.refid_parent
        
        if len(self.refid_child)>0:
            gff_defaults['child_key_transcript'] += self.refid_child
            gtf_defaults['child_key_transcript'] += self.refid_child
        
        return config_defaults, gtf_defaults, gff_defaults
    
    def integrate_assemblies_with_reference(self, locus):
        """Iterate over the set of filtered transcripts (highest TPM first)
        and update the reference by (1) replacement iff ORF is unchanged,
        (2) isoform if no containment or if ORF changes, (3) antisense, (4) intergenic.
        Updates the list of ref_reads in-place."""
        sort_order = np.argsort([-t.weight for t in locus.transcripts])
        counts_by_gene = Counter([read.attributes['gene_id'] for read in locus.ref_reads])
        for i in sort_order:
            ref_match = sense_match = antisense_match = -1
            transcript = locus.transcripts[i]
            if transcript.get_length() <= self.minlen: # Do NOT add to reference
                continue
            
            transcript.attributes['TPM'] = round(transcript.weight,1)
            ref_match = locus.matching_ref(transcript)
            if ref_match > -1:
                # Do identical work
                ref = locus.ref_reads[ref_match]
                # One final check that transcript is not truncated relative to best match ref
                left_change = ref.span[0] - transcript.span[0]
                right_change = transcript.span[1] - ref.span[1]
                if ref.attributes['source'] == 'reference': # Best match hasn't been altered yet
                    if left_change >= -self.end_cluster and right_change >= -self.end_cluster and len(ref.ranges) <= len(transcript.ranges):
                        # left and right borders did not shrink more than end_cluster
                        # Update matching ref; DO NOT add new item
                        del transcript.attributes['gene_id']
                        del transcript.attributes['transcript_id']
                        ref.attributes.update(transcript.attributes)
                        ref.ranges = transcript.ranges
                        ref.weight = transcript.weight
                        ref.capped = transcript.capped
                        ref.s_tag = transcript.s_tag
                        ref.e_tag = transcript.e_tag
                        continue
                    else: # Too much was lost on the starts/ends
                        if transcript.capped: # Pass downstream if capped
                            ref_match = -1
                        else: # Discard if uncapped (fragment)
                            continue
            
            # Continue search as a potential isoform
            sense_match = locus.sense_ref(transcript)
            if sense_match > -1:
                # Do sense overlap work
                ref = locus.ref_reads[sense_match]
                gene_id = ref.attributes['gene_id']
                counts_by_gene[gene_id] += 1
                transcript_count = counts_by_gene[gene_id]
                transcript_id = '{}.B{}'.format(gene_id, transcript_count)
                transcript.attributes['gene_id'] = gene_id
                transcript.attributes['transcript_id'] = transcript_id
            else:
                # Continue search as a potential antisense RNA
                antisense_match = locus.antisense_ref(transcript)
                if antisense_match > -1:
                    # Do antisense work
                    ref = locus.ref_reads[antisense_match]
                    gene_id = ref.attributes['gene_id']+'_AS'
                    counts_by_gene[gene_id] += 1
                    transcript_count = counts_by_gene[gene_id]
                    transcript_id = '{}.B{}'.format(gene_id, transcript_count)
                    transcript.attributes['gene_id'] = gene_id
                    transcript.attributes['transcript_id'] = transcript_id
            
            if ref_match == -1 and sense_match == -1 and antisense_match == -1:
                # Do novel transcript work
                self.new_gene_counter += 1
                gene_id = 'BOOKEND_{}'.format(self.new_gene_counter)
                counts_by_gene[gene_id] += 1
                transcript_count = counts_by_gene[gene_id] 
                transcript_id = '{}.B{}'.format(gene_id, transcript_count)
                transcript.attributes['gene_id'] = gene_id
                transcript.attributes['transcript_id'] = transcript_id

            locus.ref_reads.append(transcript)

        return sorted(locus.ref_reads)

    def output_transcripts(self, transcript, output_type):
        """Writes the RNAseqMapping object 'transcript' to an output stream,
        formatted as output_type."""
        if output_type == 'elr':
            output_line = transcript.write_as_elr()
        elif output_type == 'bed':
            output_line = transcript.write_as_bed(self.dataset.chrom_array, self.dataset.source_array, score_column='weight', name_attr='transcript_id')
        elif output_type == 'gtf':
            output_line = transcript.write_as_gtf(self.dataset.chrom_array, self.dataset.source_array[transcript.source])
        
        self.output_file.write(output_line+'\n')

    def process_entry(self, chunk):
        self.locus_counter += 1
        locus = self.make_annotation_locus(chunk)
        if locus: # Work needs to be done on non-reference transcripts
            locus.filter_fused_and_truncated_annotations()
            merged_annotations = self.integrate_assemblies_with_reference(locus)
        elif locus is not None:
            merged_annotations = locus.ref_reads
        else:
            return
        
        for transcript in merged_annotations:
            self.transcript_counter += 1
            if transcript.attributes['source'] == 'reference':
                transcript.attributes['source'] = self.refname
                transcript.source = 0
                del transcript.attributes['TPM']
            else:
                self.updated_transcript_counter += 1
                transcript.source = 1
            
            if self.fasta_out or self.orf_out:
                fasta_sequence = self.dataset.get_transcript_fasta(transcript)
                transcript.attributes['transcript_length'] = len(fasta_sequence)
            
            if self.fasta_out:
                self.output_file_fasta.write('>{}\n{}\n'.format(transcript.attributes['transcript_id'], fasta_sequence))
            
            if self.orf_out:
                orf, pos, stopless = fu.longest_orf(fasta_sequence)

                if stopless:
                    transcript.attributes['orf_length'] = 0
                else:
                    transcript.attributes['orf_length'] = len(orf)
                    self.output_file_orf.write('>{}\n{}\n'.format(transcript.attributes['transcript_id'], orf))
            
            self.output_transcripts(transcript, self.output_type)

    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend merge |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input files:\n\t{}\n".format('\n\t'.join(self.input))
        options_string += "  Reference file (-r):\n\t{}\n".format(self.reference)
        options_string += "  Output file (-o):\n\t{}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Cluster distance for ends (--end_cluster):   {}\n".format(self.end_cluster)
        options_string += "  *** Filters ***\n"
        options_string += "  Copies needed to keep assembly (--min_reps): {}\n".format(self.min_reps)
        options_string += "  Min % 5'G reads (--cap_percent):             {}\n".format(self.cap_percent)
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
        for chunk in self.generator:
            self.process_entry(chunk)
        
        self.output_file.close()
        print(self.display_summary())


if __name__ == '__main__':
    from argument_parsers import merge_parser as parser
    args = vars(parser.parse_args())
    obj = AnnotationMerger(args)
    sys.exit(obj.run())

