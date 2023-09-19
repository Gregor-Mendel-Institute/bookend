#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import copy
import bookend.core.cython_utils._rnaseq_utils as ru
from bookend.core.gtf_classify import AssemblyClassifier
from numpy import argsort
from collections import Counter, namedtuple
from math import ceil
if __name__ == '__main__':
    sys.path.append('../../bookend')

class AnnotationMerger(AssemblyClassifier):
    '''Inherits behaviors from bookend classify, but
    integrates one or more assemblies into a reference
    based on the classifications.'''
    def __init__(self, args):
        """Parses input arguments for assembly"""
        # print(args)
        self.match_data = namedtuple('match_data', 'matchtype transcript gene exonoverlap reflen tlen ref diff5p diff3p')
        self.args = args
        self.reference = self.args['REFERENCE']
        self.input = self.args['INPUT']
        self.output = self.args['OUT']
        self.end_buffer = self.args['END_BUFFER']
        self.allow_unstranded = self.args['UNSTRANDED']
        self.refname = self.args['REFNAME']
        self.min_len = self.args['MIN_LEN']
        self.min_reps = self.args['REP_FILTER']
        self.tpm_filter = self.args['TPM_FILTER']
        self.confidence = self.args['CONFIDENCE']
        self.cap_percent = self.args['CAP_PERCENT']
        self.discard = args['DISCARD']
        self.verbose = self.args['VERBOSE']
        self.attr_merge = self.args['ATTR_MERGE']
        self.keep_refs = self.args['KEEP_REFS']
        self.table = self.args['TABLE']
        self.match_unstranded = False
        if self.refname is None and self.reference is not None:
            self.refname = '.'.join(self.reference.split('/')[-1].split('.')[:-1])
        
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
            self.output_type = 'gtf'
            self.output = self.output + '.gtf'
            self.output_file = open(self.output,'w')
        
        self.table_file = open(self.table, 'w')
        self.table_file.write('transcript_id\tclass\tgene_id\tchrom\tstart\tend\tstrand\texons\tlength\tref_length\toverlap\tdiff5p\tdiff3p\tcov\tTPM\tS.reads\tE.reads\n')
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
        self.dataset.source_array = [self.refname, 'assembly']
        self.generator = self.dataset.generator
        self.locus_counter = 0
        self.new_gene_counter = 0
        self.ref_transcript_count = 0
        self.input_transcript_count = 0
        self.transcript_counter = 0
        self.updated_transcript_counter = 0
        self.removed_counter = 0
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
        self.class_colors = {
            'full_match':'29,66,134',
            'exon_match':'94,143,237',
            'isoform':'175,155,247',
            'fusion':'163,106,50',
            'fragment':'232,132,34',
            'antisense':'249,113,113',
            'intergenic':'232,187,56',
            'intronic':'174,204,106',
            'ambiguous':'200,200,200',
            'reference':'150,150,150'
        }
       
    def integrate_assemblies_with_reference(self, merged_assemblies, ref_transcripts):
        """Iterate over the set of filtered transcripts (highest TPM first)
        and update the reference by (1) replacement iff ORF is unchanged,
        (2) isoform if no containment or if ORF changes, (3) antisense, (4) intergenic.
        Updates the list of ref_reads in-place."""
        sort_order = argsort([-t.get_length() for t in merged_assemblies])
        counts_by_gene = Counter([t.attributes['gene_id'] for t in ref_transcripts])
        merged_annotations = copy.copy(ref_transcripts)
        for i in sort_order:
            transcript = merged_assemblies[i]
            if transcript.get_length() <= self.min_len: # Do NOT add to reference
                continue
            
            # Create a match_data object (matchtype transcript gene exonoverlap reflen tlen ref diff5p diff3p)
            ref_match = self.calculate_match_type(transcript, ref_transcripts)
            transcript.attributes['class'] = self.match_types[ref_match.matchtype]
            transcript.attributes['gene_id'] = ref_match.gene
            transcript.attributes['overlap'] = ref_match.exonoverlap
            transcript.attributes['ref_length'] = ref_match.reflen            
            transcript.attributes['diff5p'] = ref_match.diff5p
            transcript.attributes['diff3p'] = ref_match.diff3p
            if ref_match.matchtype == 8: # full_match
                ref = [t for t in merged_annotations if t.attributes['transcript_id'] == ref_match.transcript][0]
                self.add_rep(ref, transcript)
                ref.attributes['itemRgb'] = self.class_colors['full_match']
                continue
            elif ref_match.matchtype == 7: # exon_match, merge with ref, use .b<num> if kept
                gene_id = ref_match.transcript
            elif ref_match.matchtype in [5,6]: # fragment, fusion, treat strictly, use .b<num> if kept
                gene_id = ref_match.gene
                if not self.passes_filters(transcript, stringent=True):
                    self.removed_counter += 1
                    continue
                else: # High-confidence, keep as a valid isoform
                    gene_id = '_'.join(ref_match.gene.split(','))
            elif ref_match.matchtype in [4,1]: # isoform, use .b<num> transcript id
                gene_id = ref_match.gene
            elif ref_match.matchtype == 3: # intronic, use -IT gene suffix
                gene_id = '{}-IT'.format(ref_match.gene)
            elif ref_match.matchtype == 2: # antisense, use -AS gene suffix
                gene_id = '{}-AS'.format(ref_match.gene)
            elif ref_match.matchtype == 0: # intergenic, use BOOKEND_<num> gene name
                self.new_gene_counter += 1
                gene_id = 'BOOKEND_{}'.format(self.new_gene_counter)
            
            counts_by_gene[gene_id] += 1
            transcript_count = counts_by_gene[gene_id]
            if ref_match.matchtype == 7:
                transcript_id = '{}_{}'.format(gene_id, transcript_count)
            else:
                transcript_id = '{}.i{}'.format(gene_id, transcript_count)
            
            transcript.attributes['gene_id'] = gene_id
            transcript.attributes['transcript_id'] = transcript_id
            transcript.attributes['itemRgb'] = self.class_colors[transcript.attributes['class']]
            if '-AS' in gene_id: # Handle novel antisense isoforms
                if ref_match.matchtype in [1,4,5,6,7]:
                    transcript.attributes['class'] = 'antisense'
                    transcript.attributes['itemRgb'] = self.class_colors['antisense']
            if 'BOOKEND' in gene_id: # Handle novel intergenic isoforms
                if ref_match.matchtype in [1,4,5,6,7]:
                    transcript.attributes['class'] = 'intergenic'
                    transcript.attributes['itemRgb'] = self.class_colors['intergenic']
            
            merged_annotations.append(transcript)
            if ref_match.matchtype in [0,2]:
                ref_transcripts.append(transcript)
        
        return sorted(merged_annotations)

    def passes_filters(self, transcript, stringent=False):
        """The transcript is a suspected artifact and must pass the
        high-confidence filters defined by the user."""
        if stringent:
            multiplier = self.confidence
        else:
            multiplier = 1.0
        
        if self.keep_refs and transcript.is_reference:
            return True

        if int(transcript.attributes.get('reps',0)) < self.min_reps * multiplier:
            if self.verbose: print('Removed {} (min reps)'.format(transcript.attributes['transcript_id']))
            return False
        if not self.allow_unstranded and transcript.strand == 0:
            if self.verbose: print('Removed {} (unstranded)'.format(transcript.attributes['transcript_id']))
            return False
        if transcript.get_length() < self.min_len:
            if self.verbose: print('Removed {} (min length)'.format(transcript.attributes['transcript_id']))
            return False
        if float(transcript.attributes.get('TPM',0)) < self.tpm_filter * multiplier:
            if self.verbose: print('Removed {} (TPM filter)'.format(transcript.attributes['transcript_id']))
            return False
        if self.cap_percent > 0:
            cap_percent = float(transcript.attributes.get('S.capped',0)) / float(transcript.attributes.get('S.reads',1))
            if cap_percent < self.cap_percent:
                if self.verbose: print('Removed {} (cap percent)'.format(transcript.attributes['transcript_id']))
                return False
        
        if transcript.attributes['class'] in self.discard:
            if self.verbose: print('Removed {} (discard class: {})'.format(transcript.attributes['transcript_id'], transcript.attributes['class']))
            return False
        
        return True
    
    def output_transcript(self, transcript, output_type):
        """Writes the RNAseqMapping object 'transcript' to an output stream,
        formatted as output_type."""
        if output_type == 'elr':
            output_line = transcript.write_as_elr()
        elif output_type == 'bed':
            output_line = transcript.write_as_bed(self.dataset.chrom_array, self.dataset.source_array, score_column='weight', name_attr='transcript_id')
        elif output_type == 'gtf':
            output_line = transcript.write_as_gtf(self.dataset.chrom_array, self.dataset.source_array[transcript.source])
        
        self.output_file.write(output_line+'\n')

    def merge_assemblies(self, assemblies):
        '''Given a list of assembled transcripts,
        return a consensus set of RNAseqMapping objects.'''
        merged_assemblies = []
        # Process assemblies in order of decreasing length
        sort_order = argsort([-t.get_length() for t in assemblies])
        for i in sort_order:
            transcript = assemblies[i]
            mergematch = self.calculate_match_type(transcript, merged_assemblies)
            if mergematch.matchtype == 8: # Add attributes of assembly to existing match
                original = [t for t in assemblies if t.attributes['transcript_id'] == mergematch.transcript][0]
                self.add_rep(original, transcript)
            else:
                transcript.attributes['reps'] = 1
                merged_assemblies.append(transcript)
            
        return merged_assemblies
    
    def process_entry(self, chunk):
        self.locus_counter += 1
        ref_transcripts = [read for read in chunk if read.is_reference]
        for ref in ref_transcripts:
            ref.attributes['class'] = 'reference'
        
        nonref_transcripts = [read for read in chunk if not read.is_reference]
        self.input_transcript_count += len(nonref_transcripts)
        self.ref_transcript_count += len(ref_transcripts)
        merged_assemblies = self.merge_assemblies(nonref_transcripts)
        merged_annotations = self.integrate_assemblies_with_reference(merged_assemblies, ref_transcripts)
        for transcript in merged_annotations:
            if not self.passes_filters(transcript):
                self.removed_counter += 1
                continue
            
            if transcript.attributes.get('source','') == 'reference':
                transcript.attributes['source'] = self.refname
                transcript.attributes['itemRgb'] = self.class_colors['reference']
                transcript.source = 0
                del transcript.attributes['TPM']
            else:
                transcript.attributes['source'] = transcript.attributes['source'].replace('reference',self.refname)
                self.updated_transcript_counter += 1
            
            self.transcript_counter += 1
            self.update_table(transcript)
            self.output_transcript(transcript, self.output_type)
    
    def add_rep(self, original, rep):
        '''Incorporate the information from a repeated RNAseqMapping
        object into one that already exists. Edits original in-place.'''
        if original.attributes['source'] == 'reference': # Validation of a reference transcript
            gene_id = original.attributes['gene_id']
            transcript_id = original.attributes['transcript_id']
            original.attributes.update(rep.attributes)
            original.attributes['gene_id'] = gene_id
            original.attributes['transcript_id'] = transcript_id
            original.attributes['source'] = 'reference;{}'.format(original.attributes['source'])
            original.ranges[0] = (rep.ranges[0][0], original.ranges[0][1])
            original.ranges[-1] = (original.ranges[-1][0], rep.ranges[-1][1])
            original.span = (original.ranges[0][0], original.ranges[-1][1])
        else:
            operator = sum if self.attr_merge == 'sum' else max
            original.attributes['source'] += ';{}'.format(rep.attributes['source'])
            original.attributes['reps'] = original.attributes.get('reps',1) + 1
            if 'cov' in original.attributes:
                original.attributes['cov'] = operator([float(original.attributes['cov']), float(rep.attributes.get('cov',0))])
            
            if 'bases' in original.attributes:
                original.attributes['bases'] = operator([float(original.attributes['bases']), float(rep.attributes.get('bases',0))])
            
            if 'TPM' in original.attributes:
                original.attributes['TPM'] = operator([float(original.attributes['TPM']), float(rep.attributes.get('TPM',0))])
            
            new5p = False
            if float(rep.attributes.get('S.capped', 0)) > float(original.attributes.get('S.capped', 0)):
                # More empirical support for new 5' end
                new5p = float(rep.attributes.get('S.capped', 0)) > float(original.attributes.get('S.capped',0)) 
                original.attributes['S.left'] = min([int(original.attributes['S.left']), int(rep.attributes['S.left'])])
                original.attributes['S.right'] = max([int(original.attributes['S.right']), int(rep.attributes['S.right'])])
                original.attributes['S.reads'] = operator([float(original.attributes['S.reads']), float(rep.attributes['S.reads'])])
                original.attributes['S.capped'] = operator([float(original.attributes['S.capped']), float(rep.attributes['S.capped'])])
            
            new3p = False
            if float(rep.attributes.get('E.reads', 0)) > float(original.attributes.get('E.reads', 0)):
                # More empirical support for new 3' end
                new3p = float(rep.attributes.get('E.reads', 0)) > float(original.attributes.get('E.reads',0)) 
                original.attributes['E.left'] = min([int(original.attributes['E.left']), int(rep.attributes['E.left'])])
                original.attributes['E.right'] = max([int(original.attributes['E.right']), int(rep.attributes['E.right'])])
                original.attributes['E.reads'] = operator([float(original.attributes['E.reads']), float(rep.attributes['E.reads'])])
            
            if new5p:
                if original.strand == 1:
                    original.ranges[0] = (rep.ranges[0][0], original.ranges[0][1])
                elif original.strand == -1:
                    original.ranges[-1] = (original.ranges[-1][0], rep.ranges[-1][1])
            
            if new3p:
                if original.strand == 1:
                    original.ranges[-1] = (original.ranges[-1][0], rep.ranges[-1][1])
                elif original.strand == -1:
                    original.ranges[0] = (rep.ranges[0][0], original.ranges[0][1])
            
            original.span = (original.ranges[0][0], original.ranges[-1][1])
    
    def update_table(self, transcript):
        '''Writes a summary table like the output of bookend classify.''' 
        classification = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            transcript.attributes['transcript_id'],
            transcript.attributes['class'],
            transcript.attributes['gene_id'],
            self.dataset.chrom_array[transcript.chrom],
            transcript.ranges[0][0],
            transcript.ranges[-1][1],
            ['.','+','-'][transcript.strand],
            len(transcript.ranges),
            transcript.attributes.get('length',transcript.get_length()),
            transcript.attributes.get('ref_length','NA'),
            transcript.attributes.get('overlap','NA'),
            transcript.attributes.get('diff5p','NA'),
            transcript.attributes.get('diff3p','NA'),
            round(float(transcript.attributes.get('cov', 0)),1),
            round(float(transcript.attributes.get('TPM', 0)),1),
            round(float(transcript.attributes.get('S.reads', 0)),1),
            round(float(transcript.attributes.get('E.reads', 0)),1)
        )
        self.table_file.write(classification)

    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend merge |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input files:\n\t{}\n".format('\n\t'.join(self.input))
        options_string += "  Reference file (-r):\n\t{}\n".format(self.reference)
        options_string += "  Output file (-o):\n\t{}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Cluster distance for ends (--end_buffer):   {}\n".format(self.end_buffer)
        options_string += "  *** Filters ***\n"
        options_string += "  Discard assemblies shorter than (--min_len):   {}\n".format(self.min_len)
        options_string += "  Copies needed to keep assembly (--rep_filter): {}\n".format(self.min_reps)
        options_string += "  Discard assemblies < this TPM (--tpm_filter):  {}\n".format(self.tpm_filter)
        options_string += "  High confidence multiplier (--high_conf):      {}\n".format(self.confidence)
        options_string += "  Discard these classes (--discard):             {}\n".format(self.discard)
        options_string += "  Min % 5'G reads (--cap_percent):               {}\n".format(self.cap_percent)
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "{} loci processed ({} input transcripts, {} reference transcripts).\n".format(self.locus_counter, self.input_transcript_count, self.ref_transcript_count)
        summary += "{} transcripts removed due to filters.\n".format(self.removed_counter)
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
            if len(chunk) > 0:
                self.process_entry(chunk)
        
        self.output_file.close()
        print(self.display_summary())


if __name__ == '__main__':
    from argument_parsers import merge_parser as parser
    args = vars(parser.parse_args())
    obj = AnnotationMerger(args)
    sys.exit(obj.run())

