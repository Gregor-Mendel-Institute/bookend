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
        print(args)
        self.match_data = namedtuple('match_data', 'matchtype transcript gene exonoverlap reflen tlen ref diff5p diff3p')
        self.args = args
        self.reference = self.args['REFERENCE']
        self.input = self.args['INPUT']
        self.output = self.args['OUT']
        self.genome = self.args['GENOME']
        self.end_buffer = self.args['END_BUFFER']
        self.min_terminal = self.args['MIN_TERMINAL']
        self.allow_unstranded = self.args['UNSTRANDED']
        self.refname = self.args['REFNAME']
        self.min_len = self.args['MIN_LEN']
        self.min_start = self.args['MIN_S']
        self.min_end = self.args['MIN_E']
        self.min_reps = self.args['REP_FILTER']
        self.tpm_filter = self.args['TPM_FILTER']
        self.confidence = self.args['CONFIDENCE']
        self.cap_percent = self.args['CAP_PERCENT']
        self.discard = args['DISCARD']
        self.verbose = self.args['VERBOSE']
        self.attr_merge = self.args['ATTR_MERGE']
        self.keep_refs = self.args['KEEP_REFS']
        self.keep_best = self.args['KEEP_BEST']
        self.fusion_delim = self.args['FUSION_DELIM']
        self.table = self.args['TABLE']
        self.match_unstranded = False
        if self.refname is None and self.reference is not None:
            self.refname = '.'.join(self.reference.split('/')[-1].split('.')[:-1])
        elif self.refname is None:
            self.refname = ''
        
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
        self.table_file.write('transcript_id\tclass\tgene_id\tchrom\tstart\tend\tstrand\texons\tlength\tref_length\toverlap\tdiff5p\tdiff3p\tcov\tTPM\tS.reads\tS.capped\tE.reads\n')
        print(self.display_options())
        if self.genome is None:
            print("WARNING: no genome FASTA provided, so splice junction motifs cannot be assessed.")
        
        config_defaults, gtf_defaults, gff_defaults = self.make_config_dicts()
        self.dataset = ru.AnnotationDataset(
            annotation_files=self.input, 
            reference=self.reference,
            genome_fasta=self.genome, 
            config=config_defaults, 
            gtf_config=gtf_defaults, 
            gff_config=gff_defaults
        )
        # self.dataset.source_array = [self.refname, 'assembly']
        self.generator = self.dataset.generator
        self.locus_counter = 0
        self.new_gene_counter = 0
        self.ref_transcript_count = 0
        self.input_transcript_count = 0
        self.transcript_counter = 0
        self.updated_transcript_counter = 0
        self.removed_counter = 0
        self.match_types = {
            0:'intergenic', # 0 (lowest classification) no ref match
            1:'ambiguous',  # 1 overlapping with no strand information
            2:'antisense',  # 2 only overlaps a ref in antisense
            3:'intronic',   # 3 fully contained in a ref intron (sense)
            4:'splice_variant',    # 4 overlaps, incompatible exon chain
            4.1:'alt_tss',    # 4.1 overlaps, >=1 exon difference at 5' end
            4.2:'alt_pas',    # 4.2 overlaps, >=1 exon difference at 3' end
            4.3:'retained_intron',    # 4.3 overlaps, subset of introns
            5:'fragment',   # 5 compatible with, but fewer exons than, a ref
            6:'fusion',     # 6 shares exons with 2 or more ref genes
            7:'exon_match', # 7 shares entire exon chain, but not ends
            8:'full_match'  # 8 shares entire exon chan and ends
        }

        self.class_colors = {
            'full_match':'29,66,134',
            'exon_match':'94,143,237',
            'retained_intron':'171,182,235',
            'splice_variant':'171,103,245',
            'alt_tss':'66,157,51',
            'alt_pas':'212,38,38',
            'fusion':'137,78,23',
            'fragment':'252,116,4',
            'antisense':'255,143,143',
            'intergenic':'251,213,21',
            'intronic':'174,204,106',
            'ambiguous':'200,200,200',
            'reference':'150,150,150'
        }
       
    def integrate_assemblies_with_reference(self, merged_assemblies, ref_transcripts):
        """Iterate over the set of filtered transcripts (highest TPM first)
        and update the reference by (1) replacement iff ORF is unchanged,
        (2) isoform if no containment or if ORF changes, (3) antisense, (4) intergenic.
        Updates the list of ref_reads in-place."""
        sort_order = self.make_sort_order(merged_assemblies, type='genomic')
        counts_by_gene = Counter([t.attributes['gene_id'] for t in ref_transcripts])
        merged_annotations = copy.copy(ref_transcripts)
        for i in sort_order:
            transcript = merged_assemblies[i]
            # Create a match_data object (matchtype transcript gene exonoverlap reflen tlen ref diff5p diff3p)
            ref_match = self.calculate_match_type(transcript, merged_annotations)
            transcript.attributes['class'] = self.match_types[ref_match.matchtype]
            transcript.attributes['gene_id'] = ref_match.gene
            transcript.attributes['overlap'] = ref_match.exonoverlap
            transcript.attributes['ref_length'] = ref_match.reflen
            transcript.attributes['diff5p'] = ref_match.diff5p
            transcript.attributes['diff3p'] = ref_match.diff3p
            if ref_match.matchtype >= 4:
                try:
                    ref = [t for t in merged_annotations if t.attributes['transcript_id'] == ref_match.transcript][0]
                    # ref.attributes['traceback'] = ref.attributes['transcript_id']
                except:
                    ref = None
            
            if ref_match.matchtype == 8: # full_match
                self.add_rep(ref, transcript)
                ref.attributes['itemRgb'] = self.class_colors['full_match']
                continue
            elif ref_match.matchtype == 7: # exon_match, merge with ref, use .i<num> if kept
                self.add_rep(ref, transcript)
                continue
            elif ref_match.matchtype in [5,6]: # fragment, fusion, treat strictly, use .i<num> if kept
                gene_id = ref_match.gene
            elif int(ref_match.matchtype) in [4,1]: # isoform, use .i<num> transcript id
                gene_id = ref_match.gene
            elif ref_match.matchtype == 3: # intronic, use -IT gene suffix
                gene_id = '{}-IT'.format(ref_match.gene)
            elif ref_match.matchtype == 2: # antisense, use -AS gene suffix
                gene_id = '{}-AS'.format(ref_match.gene)
            elif ref_match.matchtype == 0: # intergenic, use BOOKEND_<num> gene name
                self.new_gene_counter += 1
                gene_id = 'BOOKEND-{}'.format(self.new_gene_counter)
            
            
            gene_id = self.fusion_delim.join(sorted(set(gene_id.split(self.fusion_delim))))
            
            counts_by_gene[gene_id] += 1
            transcript_count = counts_by_gene[gene_id]
            transcript_id = '{}.i{}'.format(gene_id, transcript_count)
            
            transcript.attributes['gene_id'] = gene_id
            transcript.attributes['transcript_id'] = transcript_id
            transcript.attributes['itemRgb'] = self.class_colors[transcript.attributes['class']]
            
            if 4 <= ref_match.matchtype < 8 and ref is not None:
                operator = sum if self.attr_merge == 'sum' else max
                self.merge_ends(ref, transcript, operator, both=True)
                new_match = self.calculate_match_type(transcript, [ref])
                if new_match.matchtype == 8:
                    self.add_rep(ref, transcript)
                else:
                    merged_annotations.append(transcript)
            else:
                merged_annotations.append(transcript)
        
        return sorted(merged_annotations)
    
    def noncanonical_junctions(self, transcript):
        canonical = ''
        if transcript.strand == 1:
            canonical = ['GTAG','GCAG']
        elif transcript.strand == -1:
            canonical = ['CTAC', 'CTGC']
        
        for left,right in transcript.junctions():
            motif = (
                ru.get_flank(self.dataset.genome, self.dataset.chrom_array[transcript.chrom], left-1, 1, 'E', 2) + 
                ru.get_flank(self.dataset.genome, self.dataset.chrom_array[transcript.chrom], right, 1, 'S', 2)
            ).upper()
            if motif not in canonical:
                return True
        
        return False
    
    def passes_filters(self, transcript):
        """The transcript is a suspected artifact and must pass the
        high-confidence filters defined by the user."""
        robust = transcript.attributes.get('robust',False)
        stringent = transcript.attributes.get('class', 'reference') in ['fusion','fragment','retained_intron']
        if self.dataset.has_genome:
            stringent = stringent or self.noncanonical_junctions(transcript)
        
        if stringent:
            multiplier = self.confidence
        else:
            multiplier = 1.0

        if self.keep_refs and transcript.is_reference:
            return True

        min_reps = 1 if robust else self.min_reps
        if transcript.attributes['class'] in self.discard:
            if self.verbose: print('Removed {} (discard class: {})'.format(transcript.attributes['transcript_id'], transcript.attributes['class']))
            return False

        if int(transcript.attributes.get('reps',0)) < min_reps * multiplier:
            if self.verbose: print('Removed {} (min reps)'.format(transcript.attributes['transcript_id']))
            return False
        if not self.allow_unstranded and transcript.strand == 0:
            if self.verbose: print('Removed {} (unstranded)'.format(transcript.attributes['transcript_id']))
            return False
        if transcript.get_length() < self.min_len:
            if self.verbose: print('Removed {} (min length)'.format(transcript.attributes['transcript_id']))
            return False
        if (transcript.ranges[0][1]-transcript.ranges[0][0]) < self.min_terminal or (transcript.ranges[-1][1]-transcript.ranges[-1][0]) < self.min_terminal:
            if int(transcript.attributes.get('reps',0)) < min_reps * self.confidence:
                if self.verbose: print('Removed {} (min terminal exon)'.format(transcript.attributes['transcript_id']))
                return False
        if float(transcript.attributes.get('TPM',0)) < self.tpm_filter * multiplier:
            if self.verbose: print('Removed {} (TPM filter)'.format(transcript.attributes['transcript_id']))
            return False
        if transcript.attributes.get('class', 'reference') in ['fragment','alt_tss']:
            cap_percent = float(transcript.attributes.get('S.capped',0)) / float(transcript.attributes.get('S.reads',1))
            if cap_percent < self.cap_percent:
                if self.verbose: print('Removed {} (cap percent)'.format(transcript.attributes['transcript_id']))
                return False

            if self.min_start * self.confidence > float(transcript.attributes.get('S.reads',self.min_start)):
                if self.verbose: print('Removed {} (min start reads)'.format(transcript.attributes['transcript_id']))
                return False

        if self.min_start > float(transcript.attributes.get('S.reads',self.min_start)):
            if self.verbose: print('Removed {} (min start reads)'.format(transcript.attributes['transcript_id']))
            return False

        if transcript.attributes.get('class', 'reference') == 'alt_pas':
            if self.min_end * self.confidence > float(transcript.attributes.get('E.reads',self.min_end)):
                if self.verbose: print('Removed {} (min end reads)'.format(transcript.attributes['transcript_id']))
                return False
        else:
            if self.min_end > float(transcript.attributes.get('E.reads',self.min_end)):
                if self.verbose: print('Removed {} (min end reads)'.format(transcript.attributes['transcript_id']))
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
        elif 'gff' in output_type:
            output_line = transcript.write_as_gtf(self.dataset.chrom_array, self.dataset.source_array[transcript.source],attr_format='GFF')
        
        self.output_file.write(output_line+'\n')

    def merge_assemblies(self, assemblies):
        '''Given a list of assembled transcripts,
        return a consensus set of RNAseqMapping objects.'''
        merged_assemblies = []
        sort_order = self.make_sort_order(assemblies, type='genomic')
        for i in sort_order:
            transcript = assemblies[i]
            samples = int(transcript.attributes.get('samples',-1))
            if samples == 0:
                if self.verbose:
                    if self.verbose: print('Removed {} (not expressed)'.format(transcript.attributes['transcript_id']))
                
                continue
            
            transcript.attributes['traceback'] = transcript.attributes['transcript_id']
            transcript.attributes['robust'] = True if samples >= 3 else False
            mergematch = self.calculate_match_type(transcript, merged_assemblies)
            if mergematch.matchtype == 8: # Add attributes of assembly to existing match
                original = [t for t in assemblies if t.attributes['transcript_id'] == mergematch.transcript][0]
                self.add_rep(original, transcript)
            else:
                transcript.attributes['reps'] = 1
                merged_assemblies.append(transcript)
            
        return merged_assemblies
    
    def exon_count(self, transcript, min_size=1):
        """Returns the number of exons in a transcript."""
        return sum([r-l >= min_size for l,r in transcript.ranges])

    def make_sort_order(self, transcripts, type='exonic'):
        """Returns a sort order for the given list of transcripts
        based on exonic length weighted by support."""
        if type == 'exonic':
            sort_order = argsort([-(float(t.get_length()) * float(t.weight) * float(self.exon_count(t, 20))) for t in transcripts])
        else: # genomic
            sort_order = argsort([-((t.ranges[-1][1] - t.ranges[0][0]) * float(t.weight) * float(self.exon_count(t, 20))) for t in transcripts])
        
        return sort_order

    def process_entry(self, chunk):
        self.locus_counter += 1
        for transcript in chunk:
            transcript.attributes['locus_size'] = len(chunk)
        
        ref_transcripts = [read for read in chunk if read.is_reference]
        for ref in ref_transcripts:
            ref.attributes['class'] = 'reference'
        
        nonref_transcripts = [read for read in chunk if not read.is_reference]
        self.input_transcript_count += len(nonref_transcripts)
        self.ref_transcript_count += len(ref_transcripts)
        merged_assemblies = self.merge_assemblies(nonref_transcripts)
        merged_annotations = self.integrate_assemblies_with_reference(merged_assemblies, ref_transcripts)
        sort_order = self.make_sort_order(merged_annotations, type='genomic')
        for i in sort_order:
            transcript = merged_annotations[i]
            if self.keep_best and i == sort_order[0]:
                pass
            else:
                if not self.passes_filters(transcript):
                    self.removed_counter += 1
                    continue                
            
            if transcript.attributes.get('soure','') == 'reference':
                transcript.attributes['source'] = self.refname
                transcript.attributes['itemRgb'] = self.class_colors['reference']
                transcript.source = 0
                del transcript.attributes['TPM']
            elif transcript.attributes.get('source','') == '':
                transcript.attributes['source'] = self.dataset.source_array[transcript.source]
            else:
                transcript.attributes['source'] = transcript.attributes['source'].replace('reference',self.refname)
                self.updated_transcript_counter += 1
            
            if 'S.peak' in transcript.attributes:
                del transcript.attributes['S.peak']
            
            if 'E.peak' in transcript.attributes:
                del transcript.attributes['E.peak']
            
            transcript.attributes['S.reads'] = round(float(transcript.attributes.get('S.reads',0)),1)
            transcript.attributes['S.capped'] = round(float(transcript.attributes.get('S.capped',0)),1)
            transcript.attributes['E.reads'] = round(float(transcript.attributes.get('E.reads',0)),1)
            transcript.attributes['bases'] = round(float(transcript.attributes.get('bases',0)))
            transcript.attributes['TPM'] = round(float(transcript.attributes.get('TPM',0)),2)
            self.transcript_counter += 1
            self.update_table(transcript)
            self.output_transcript(transcript, self.output_type)
    
    def merge_ends(self, t1, t2, operator, both=False):
        """Combine end information of two compatible transcripts.
        Modifies t1 (and t2 if 'both') in-place to include the 
        full TSS and PAS ranges."""
        if t1.strand != t2.strand:
            return
        
        replace = t1.is_reference and not both # This transcript had no end support, ignore it
        t1_sl = int(t1.attributes.get('S.left', -self.end_buffer))
        t1_sr = int(t1.attributes.get('S.right', -self.end_buffer))
        t1_el = int(t1.attributes.get('E.left', -self.end_buffer))
        t1_er = int(t1.attributes.get('E.right', -self.end_buffer))
        t2_sl = int(t2.attributes.get('S.left', self.end_buffer))
        t2_sr = int(t2.attributes.get('S.right', self.end_buffer))
        t2_el = int(t2.attributes.get('E.left', self.end_buffer))
        t2_er = int(t2.attributes.get('E.right', self.end_buffer))
        if t1.strand == 1:
            t1_spos = t1.ranges[0][0]
            t1_epos = t1.ranges[-1][1]
            t2_spos = t2.ranges[0][0]
            t2_epos = t2.ranges[-1][1]
        elif t1.strand == -1:
            t1_spos = t1.ranges[-1][1]
            t1_epos = t1.ranges[0][0]
            t2_spos = t2.ranges[-1][1]
            t2_epos = t2.ranges[0][0]
        
        # Check that ranges overlap
        new5p = replace or (t1_sl <= t2_sr and t1_sr >= t2_sl) or self.end_buffer >= abs(t1_spos - t2_spos)
        new3p = replace or (t1_el <= t2_er and t1_er >= t2_el) or self.end_buffer >= abs(t1_epos - t2_epos)
        # Check if exon boundaries prevent merge
        if t1.strand == 1:
            new5p = new5p and (t2_sr <= (t1.ranges[0][1] - self.min_terminal) and t1_sr <= (t2.ranges[0][1] - self.min_terminal))
            new3p = new3p and (t2_el >= (t1.ranges[-1][0] + self.min_terminal) and t1_el >= (t2.ranges[-1][0] + self.min_terminal))
        elif t1.strand == -1:
            new5p = new5p and (t2_sl >= (t1.ranges[-1][0] + self.min_terminal) and t1_sl >= (t2.ranges[-1][0] + self.min_terminal))
            new3p = new3p and (t2_er <= (t1.ranges[0][1] - self.min_terminal) and t1_er <= (t2.ranges[0][1] - self.min_terminal))
        
        if new5p: # The 5' end needs to be updated
            t1_s = float(t1.attributes.get('S.peak',t1.attributes.get('S.reads',0)))
            t2_s = float(t2.attributes.get('S.peak',t2.attributes.get('S.reads',0)))
            if t1_s > t2_s:
                spos = t1_spos
                t1.attributes['S.peak'] = t1_s
                t2.attributes['S.peak'] = t1_s
                s_left = t2_sl if replace else min([t1_sl, t2_spos])
                s_right = t2_sr if replace else max([t1_sr, t2_spos])
            else:
                spos = t2_spos
                t1.attributes['S.peak'] = t2_s
                t2.attributes['S.peak'] = t2_s
                s_left = t2_sl if replace else min([t2_sl, t1_spos])
                s_right = t2_sr if replace else max([t2_sr, t1_spos])
            
            t1.attributes['S.left'] = s_left
            t1.attributes['S.right'] = s_right
            if both:
                t2.attributes['S.left'] = s_left
                t2.attributes['S.right'] = s_right
            
            if t1.strand == 1:
                t1.ranges[0] = (spos, t1.ranges[0][1])
                if both:
                    t2.ranges[0] = (spos, t2.ranges[0][1])
            elif t1.strand == -1:
                t1.ranges[-1] = (t1.ranges[-1][0], spos)
                if both:
                    t2.ranges[-1] = (t2.ranges[-1][0], spos)
        
        if new3p: # The 3' end nees to be updated
            t1_e = float(t1.attributes.get('E.peak',t1.attributes.get('E.reads',0)))
            t2_e = float(t2.attributes.get('E.peak',t2.attributes.get('E.reads',0)))
            if t1_e > t2_e:
                epos = t1_epos
                t1.attributes['E.peak'] = t1_e
                t2.attributes['E.peak'] = t1_e
                e_left = t1_el if replace else min([t1_el, t2_epos])
                e_right = t1_er if replace else max([t1_er, t2_epos])
            else:
                epos = t2_epos
                t1.attributes['E.peak'] = t2_e
                t2.attributes['E.peak'] = t2_e
                e_left = t2_el if replace else min([t2_el, t1_epos])
                e_right = t2_er if replace else max([t2_er, t1_epos])
            
            t1.attributes['E.left'] = e_left
            t1.attributes['E.right'] = e_right
            if both:
                t2.attributes['E.left'] = e_left
                t2.attributes['E.right'] = e_right
            
            if t1.strand == 1:
                t1.ranges[-1] = (t1.ranges[-1][0], epos)
                if both:
                    t2.ranges[-1] = (t2.ranges[-1][0], epos)
            elif t1.strand == -1:
                t1.ranges[0] = (epos, t1.ranges[0][1])
                if both:
                    t2.ranges[0] = (epos, t2.ranges[0][1])
        
        t1.span = (t1.ranges[0][0], t1.ranges[-1][1])
        t2.span = (t2.ranges[0][0], t2.ranges[-1][1])


    def add_rep(self, original, rep):
        '''Incorporate the information from a repeated RNAseqMapping
        object into one that already exists. Edits original in-place.'''
        operator = sum if self.attr_merge == 'sum' else max
        if original.attributes.get('class','') == 'reference': # Validation of a reference transcript
            gene_id = original.attributes['gene_id']
            transcript_id = original.attributes['transcript_id']
            original.attributes.update(rep.attributes)
            original.attributes['gene_id'] = gene_id
            original.attributes['transcript_id'] = transcript_id
            original.attributes['source'] = 'reference;{}'.format(original.attributes['source'])
            original.attributes['robust'] = rep.attributes.get('robust',False)
            if 'traceback' in original.attributes:
                original.attributes['traceback'] += ';'+rep.attributes['transcript_id']
            else:
                original.attributes['traceback'] = original.attributes['transcript_id']+';'+rep.attributes['transcript_id']
            
            original.ranges[0] = (rep.ranges[0][0], original.ranges[0][1])
            original.ranges[-1] = (original.ranges[-1][0], rep.ranges[-1][1])
            original.span = (original.ranges[0][0], original.ranges[-1][1])
            if 'samples' in rep.attributes: original.attributes['samples'] = operator([int(rep.attributes['samples']), int(rep.attributes.get('samples',1))])
            if 'cov' in rep.attributes: original.attributes['cov'] = operator([float(rep.attributes['cov']), float(rep.attributes.get('cov',0))])
            if 'bases' in rep.attributes: original.attributes['bases'] = operator([float(rep.attributes['bases']), float(rep.attributes.get('bases',0))])
            if 'TPM' in rep.attributes: original.attributes['TPM'] = operator([float(rep.attributes['TPM']), float(rep.attributes.get('TPM',0))])
            if 'S.reads' in rep.attributes: original.attributes['S.reads'] = operator([float(rep.attributes['S.reads']), float(rep.attributes.get('S.reads',0))])
            if 'S.capped' in rep.attributes: original.attributes['S.capped'] = operator([float(rep.attributes['S.capped']), float(rep.attributes.get('S.capped',0))])
            if 'E.reads' in rep.attributes: original.attributes['E.reads'] = operator([float(rep.attributes['E.reads']), float(rep.attributes.get('E.reads',0))])
            if 'S.left' in rep.attributes: original.attributes['S.left'] = int(rep.attributes['S.left'])
            if 'S.right' in rep.attributes: original.attributes['S.right'] = int(rep.attributes['S.right'])
            if 'E.left' in rep.attributes: original.attributes['E.left'] = int(rep.attributes['E.left'])
            if 'E.right' in rep.attributes: original.attributes['E.right'] = int(rep.attributes['E.right'])            
        else:
            original.attributes['robust'] = original.attributes.get('robust',False) or rep.attributes.get('robust',False)
            if 'samples' in original.attributes: original.attributes['samples'] = operator([int(original.attributes['samples']), int(rep.attributes.get('samples',1))])
            if 'cov' in original.attributes: original.attributes['cov'] = operator([float(original.attributes['cov']), float(rep.attributes.get('cov',0))])
            if 'bases' in original.attributes: original.attributes['bases'] = operator([float(original.attributes['bases']), float(rep.attributes.get('bases',0))])
            if 'TPM' in original.attributes: original.attributes['TPM'] = operator([float(original.attributes['TPM']), float(rep.attributes.get('TPM',0))])
            if 'S.reads' in original.attributes: original.attributes['S.reads'] = operator([float(original.attributes['S.reads']), float(rep.attributes.get('S.reads',0))])
            if 'S.capped' in original.attributes: original.attributes['S.capped'] = operator([float(original.attributes['S.capped']), float(rep.attributes.get('S.capped',0))])
            if 'E.reads' in original.attributes: original.attributes['E.reads'] = operator([float(original.attributes['E.reads']), float(rep.attributes.get('E.reads',0))])            
            self.merge_ends(original, rep, operator)
        
        reps = sorted(set(original.attributes['source'].split(';')+rep.attributes['source'].split(';')))
        if 'traceback' in original.attributes:
            original.attributes['traceback'] += ';'+rep.attributes['transcript_id']
        else:
            original.attributes['traceback'] = original.attributes['transcript_id']+';'+rep.attributes['transcript_id']
        
        original.attributes['source']  =  ';'.join(reps)
        original.attributes['reps'] = len(reps)
    
    def update_table(self, transcript):
        '''Writes a summary table like the output of bookend classify.''' 
        classification = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
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
            round(float(transcript.attributes.get('S.capped', 0)),1),
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

