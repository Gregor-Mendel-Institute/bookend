#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import time
import gzip
import pysam
from multiprocessing.dummy import Pool as ThreadPool    
from bookend.core.cython_utils._rnaseq_utils import RNAseqDataset, read_generator, config_defaults
from bookend.core.cython_utils._assembly_utils import Locus
from bookend.core.elr_combine import ELRcombiner

if __name__ == '__main__':
    sys.path.append('../../bookend')


class Assembler:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        print(args)
        self.start_time = time.time()
        self.output = args['OUT']
        self.source = args['SOURCE']
        self.cov_out = args['COV_OUT']
        self.incomplete = args['INCOMPLETE']
        self.max_gap = args['MAX_GAP']
        self.end_cluster = args['END_CLUSTER']
        self.min_overhang = args['MIN_OVERHANG']
        self.min_cov = args['MIN_COV']
        self.min_unstranded_cov = args['MIN_UNSTRANDED']
        self.min_start = args['MIN_S']
        self.min_end = args['MIN_E']
        self.min_intron_length = args['MIN_INTRON_LEN']
        self.intron_filter = args['INTRON_FILTER']
        self.truncation_filter = args['TRUNCATION_FILTER']
        self.min_proportion = args['MIN_PROPORTION']
        self.cap_bonus = args['CAP_BONUS']
        self.cap_filter = args['CAP_FILTER']
        self.minlen = args['MINLEN']
        self.verbose = args['VERBOSE']
        self.input = args['INPUT']
        self.ignore_labels = args['IGNORE_LABELS']
        self.ignore_sources = not args['USE_SOURCES']
        self.require_cap = args['REQUIRE_CAP']
        self.antisense_filter = 0.01
        config = config_defaults
        if self.ignore_labels:
            self.incomplete = True
        
        if self.incomplete: # Some settings are incompatible with writing incomplete transcripts.
            self.min_start = 0
            self.min_end = 0
        
        if len(self.input) == 1:
            self.input = self.input[0]
            if self.input_is_valid(self.input):
                self.file_type = self.file_extension(self.input)
                if self.file_type in ['bam','sam']:
                    save = pysam.set_verbosity(0)
                    self.input_file = pysam.AlignmentFile(self.input)
                    save = pysam.set_verbosity(save)
                    try:
                        longread = 'minimap' in self.input_file.header['PG'][0]['ID']
                    except:
                        longread = False
                    
                    if longread:
                        print("WARNING: Detected long-read alignments. Setting defaults for unlabeled long-read cDNA.")
                        print("Read the Bookend User Guide for long-read analysis recommendations.")
                        long_defaults = {
                            'source':'',
                            's_tag':True,
                            'e_tag':True,
                            'labels_are_trimmed':False,
                            'capped':False,
                            'stranded':False,
                            'reverse':False,
                            'start_seq':None,
                            'end_seq':None,
                            'mismatch_rate':0.2,
                            'error_rate':0.2,
                            'min_reps':1,
                            'cap_bonus':1,
                            'sj_shift':2,
                            'max_headclip':120,
                            'gene_delim':'.',
                            'remove_noncanonical':False,
                            'secondary':False,
                            'ignore_ends':False,
                            'quality_filter':True,
                            'reference':None,
                            'sj':None,
                        }
                        config.update(long_defaults)
                    
                    self.dataset = RNAseqDataset(chrom_array=self.input_file.header.references, config=config)
                elif self.file_type == 'elr.gz':
                    self.dataset = RNAseqDataset()
                    self.input_file = gzip.open(self.input, 'rt')
                else:
                    self.dataset = RNAseqDataset()
                    self.input_file = open(self.input, 'r')
            else:
                print("\nERROR: input file must be a valid format (BED, ELR, BAM, SAM).")
                sys.exit(1)
        elif len(self.input) > 1: # Interleave multiple input files for assembly
            if not all([self.input_is_valid(filename, valid_formats=['elr','elr.gz']) for filename in self.input]):
                print("\nERROR: Multi-input assembly can only be performed on position-sorted ELR files.")
                sys.exit(1)
            
            self.dataset = RNAseqDataset()
            self.file_type = 'elr'
            combine_args = {
                'INPUT':self.input,
                'OUTPUT':'stdout',
                'TEMPDIR':'{}_combinetmp'.format(self.input[0])
            }
            combiner = ELRcombiner(combine_args)
            self.input_file = combiner.combine_files(combiner.input, combiner.output_file, iterator=True)
            header = next(self.input_file)
            self.dataset.chrom_array = combiner.dataset.chrom_array
            self.dataset.source_array = combiner.dataset.source_array
            self.dataset.chrom_dict = combiner.dataset.chrom_dict
            self.dataset.source_dict = combiner.dataset.source_dict
        else:
            print("\nERROR: No input file(s) provided.")
            sys.exit(1)
        
        self.complete = not self.incomplete
        self.output_type = self.file_extension(self.output)
        if self.output_type is None:
            self.output_type = 'gtf'
        
        self.generator = read_generator(self.input_file, self.dataset, self.file_type, self.max_gap, 0)
        self.chunk_counter = 0
        self.transcript_counter = 0
        self.output_file = open(self.output,'w')
    
    def output_transcripts(self, transcript, output_type):
        """Writes the RNAseqMapping object 'transcript' to an output stream,
        formatted as output_type."""
        if output_type == 'elr':
            output_line = transcript.write_as_elr(endweights=True)
        elif output_type == 'bed':
            output_line = transcript.write_as_bed(self.dataset.chrom_array, ['bookend'], score_column='coverage')
        elif output_type == 'gtf':
            output_line = transcript.write_as_gtf(self.dataset.chrom_array, 'bookend')
        
        self.output_file.write(output_line+'\n')
    
    def process_entry(self, chunk):
        STOP_AT=float('inf')
        # STOP_AT=1000000
        if len(chunk) > 0:            
            chrom = chunk[0].chrom
            self.chunk_counter += 1
            if self.verbose:
                print("\nProcessing chunk {} ({}:{}-{}, {} reads)".format(
                    self.chunk_counter,
                    self.dataset.chrom_array[chrom], chunk[0].left(),chunk[-1].right(),
                    len(chunk)
                ), end=" ")
            
            locus = Locus(
                chrom=chrom, 
                chunk_number=self.chunk_counter, 
                list_of_reads=chunk,
                max_gap=self.max_gap,
                end_cluster=self.end_cluster,
                min_overhang=self.min_overhang, 
                reduce=True, 
                minimum_proportion=self.min_proportion, 
                min_intron_length=self.min_intron_length, 
                antisense_filter=self.antisense_filter, 
                cap_bonus=self.cap_bonus,
                cap_filter=self.cap_filter, 
                complete=False, 
                verbose=self.verbose, 
                naive=self.ignore_sources, 
                intron_filter=self.intron_filter,
                splittable=True,
                ignore_ends=self.ignore_labels, 
                allow_incomplete=self.incomplete,
                require_cap=self.require_cap,
                min_start=self.min_start,
                min_end=self.min_end,
                assemble=True,
                truncation_filter=self.truncation_filter,
            )
            self.chunk_counter = locus.chunk_number
            total_bases = locus.bases
            transcripts_written = 0
            if total_bases > 0:
                bases_used = 0
                if self.verbose:
                    print('\n[{}:{}-{}] '.format(self.dataset.chrom_array[chrom], chunk[0].left(),chunk[-1].right()), end=" ")
                
                for transcript in locus.transcripts:
                    if self.passes_all_checks(transcript):
                        self.output_transcripts(transcript, self.output_type)
                        if self.cov_out:
                            if self.ignore_sources:
                                self.covfile.write('{}\t{}\n'.format(transcript.attributes['transcript_id'], round(transcript.coverage, 1)))
                            else:
                                source_cov = [0.]*len(self.dataset.source_array)
                                for k,v in locus.assembly_source_cov[transcript.attributes['transcript_id']].items():
                                    source_cov[k] = v
                                
                                self.covfile.write('{}\t{}\n'.format(transcript.attributes['transcript_id'], '\t'.join([str(round(v,1)) for v in source_cov])))
                        
                        bases_used += transcript.attributes['bases']
                        transcripts_written += 1
                
                if self.verbose:
                    print('{} transcripts from {}/{} bases ({}%)'.format(
                        transcripts_written, round(bases_used,1), round(total_bases,1), round(bases_used/total_bases*100,2)), end=" ")
                
                if chunk[0].left() >= STOP_AT:
                    sys.exit()
                
                self.transcript_counter += transcripts_written
            
            del locus
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend assemble |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                                       {}\n".format(self.input)
        options_string += "  Output file (-o):                                 {}\n".format(self.output)
        options_string += "  Source name (--source):                           {}\n".format(self.source)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Max allowed gap in coverage (--max_gap):          {}\n".format(self.max_gap)
        options_string += "  Max end cluster distance (--end_cluster):         {}\n".format(self.end_cluster)
        options_string += "  Min spanning bases (--min_overhang):              {}\n".format(self.min_overhang)
        options_string += "  Split read sources (--use_sources):               {}\n".format(not self.ignore_sources)
        options_string += "  Ignore end labels (--ignore_labels):              {}\n".format(self.ignore_labels)
        options_string += "  *** Filters ***\n"
        options_string += "  Min bp transcript length (--min_len):             {}\n".format(self.minlen)
        options_string += "  Min isoform coverage (--min_cov):                 {}\n".format(self.min_cov)
        options_string += "  Min unstranded coverage (--min_unstranded_cov):   {}\n".format(self.min_unstranded_cov)
        options_string += "  Min isoform contribution (--min_proportion):      {}\n".format(self.min_proportion)
        options_string += "  Min retained intron proportion (--intron_filter): {}\n".format(self.intron_filter)
        options_string += "  Min number of Start Tags (--min_start):           {}\n".format(self.min_start)
        options_string += "  Min number of End Tags (--min_end):               {}\n".format(self.min_end)
        options_string += "  Min percent 5' end reads w/ uuG (--cap_filter):   {}\n".format(self.cap_filter)
        options_string += "  Score multiplier for Cap Tags (--cap_bonus):      {}\n".format(self.cap_bonus)
        options_string += "  Keep fragmented assemblies (--allow_incomplete):  {}\n".format(self.incomplete)
        return options_string
    
    def display_summary(self):
        summary = ''
        summary += '\nTotal elapsed time: {}'.format(round(self.end_time - self.start_time, 5))
        summary += '\n{} assembled transcripts written to {}'.format(self.transcript_counter, self.output)
        return summary
    
    def file_extension(self, filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        elif split_name[-1].lower() == 'gz':
            return '.'.join([n.lower() for n in split_name[-2:]])
        else:
            return split_name[-1].lower()
    
    def input_is_valid(self, filename, valid_formats=['bed','elr','bam','sam','gtf','gff3','gff','elr.gz']):
        """Boolean if the file is a format that Assembler can parse."""
        if self.file_extension(filename) in valid_formats:
            return True
        else:
            return False
    
    def passes_all_checks(self, transcript):
        """Determines if a transcript model passes the collection
        of filters set by the commandline arguments.
        'transcript' is an RNAseqMapping object,
        see _rnaseq_utils.pyx for details.
        """
        if transcript.coverage < self.min_cov: return False
        if not self.incomplete and not transcript.complete: return False
        if transcript.attributes['length'] < self.minlen: return False
        if transcript.attributes.get('S.reads',0) < self.min_start: return False
        if transcript.attributes.get('E.reads',0) < self.min_end: return False
        # if not args.INCOMPLETE and True not in transcript.splice and transcript.attributes['S.capped']/transcript.attributes['S.reads'] < args.CAP_PERCENT: return False
        if transcript.strand == 0 and transcript.coverage < self.min_unstranded_cov: return False
        return True
    
    def run(self):
        """Executes end labeling on all reads."""
        if self.output != 'stdout':
            print(self.display_options())
        
        if self.cov_out:
            self.covfile=open(self.cov_out, 'w')

        if self.output_type != 'gtf':
            self.output_file.write('\n'.join(self.dataset.dump_header())+'\n')
        
        if self.cov_out:
            if not self.ignore_sources:
                self.covfile.write('{}\n'.format('\t'.join(self.dataset.source_array)))
        
        for chunk in self.generator:
            self.process_entry(chunk)
        
        if len(self.input) == 1:
            self.output_file.close()
        
        if self.cov_out:
            self.covfile.close()
        
        self.end_time = time.time()
        print(self.display_summary())

if __name__ == '__main__':
    from argument_parsers import assemble_parser as parser
    args = vars(parser.parse_args())
    obj = Assembler(args)
    sys.exit(obj.run())

