#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import time
import pysam
import bookend.core.cython_utils._fasta_utils as fu
import bookend.core.cython_utils._rnaseq_utils as ru
import bookend.core.cython_utils._assembly_utils as au
if __name__ == '__main__':
    sys.path.append('../../bookend')


class Assembler:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.start_time = time.time()
        self.output = args['OUT']
        self.source = args['SOURCE']
        self.genome = args['GENOME']
        self.incomplete = args['INCOMPLETE']
        self.infer_starts = args['INFER_STARTS']
        self.infer_ends = args['INFER_ENDS']
        self.max_gap = args['MAX_GAP']
        self.end_cluster = args['END_CLUSTER']
        self.min_overhang = args['MIN_OVERHANG']
        self.min_cov = args['MIN_COV']
        self.min_unstranded_cov = args['MIN_UNSTRANDED']
        self.min_start = args['MIN_S']
        self.min_end = args['MIN_E']
        self.intron_filter = args['INTRON_FILTER']
        self.min_proportion = args['MIN_PROPORTION']
        self.cap_percent = args['CAP_PERCENT']
        self.minlen = args['MINLEN']
        self.verbose = args['VERBOSE']
        self.input = args['INPUT']
        
        if self.incomplete: # Some settings are incompatible with writing incomplete transcripts.
            self.infer_starts = self.infer_ends = True
        
        if self.infer_starts:
            self.min_start = 0
        
        if self.infer_ends:
            self.min_end = 0
        
        if self.input_is_valid(self.input):
            self.file_type = self.file_extension(self.input)
            if self.file_type in ['bam','sam']:
                self.input_file = pysam.AlignmentFile(self.input)
                self.dataset = ru.RNAseqDataset(genome_fasta=self.genome, chrom_array=self.input_file.header.references)
            else:
                self.dataset = ru.RNAseqDataset()
                self.input_file = open(self.input)
        else:
            print("\nERROR: input file must be a valid format (BED, ELR, BAM, SAM).")
            sys.exit(1)
        
        self.complete = not self.incomplete
        self.output_type = self.file_extension(self.output)
        if self.output_type is None:
            self.output_type = 'gtf'
        
        self.generator = ru.read_generator(self.input_file, self.dataset, self.file_type, self.max_gap)
        self.chunk_counter = 0
        self.output_file = open(self.output,'w')
    
    def output_transcripts(self, transcript, output_type):
        """Writes the RNAseqMapping object 'transcript' to an output stream,
        formatted as output_type."""
        if output_type == 'elr':
            output_line = transcript.write_as_elr()
        elif output_type == 'bed':
            output_line = transcript.write_as_bed(self.dataset.chrom_array, ['bookend'], score_column='coverage')
        elif output_type == 'gtf':
            output_line = transcript.write_as_gtf(self.dataset.chrom_array, 'bookend')
        
        self.output_file.write(output_line+'\n')

    def process_entry(self, chunk):
        STOP_AT=float('inf')
        # STOP_AT=3624
        if len(chunk) > 0:
            chrom = chunk[0].chrom
            subchunks = ru.split_chunk(chunk, self.min_proportion, self.max_gap)
            for subchunk in subchunks:
                self.chunk_counter += 1
                if len(subchunk) > 0:
                    if self.verbose:
                        print('\n[{}:{}-{}] Processing chunk.'.format(self.dataset.chrom_array[chrom], subchunk[0].left(),subchunk[-1].right()))
                    
                    locus = au.Locus(chrom, self.chunk_counter, subchunk, self.max_gap, self.end_cluster, self.min_overhang, True, self.min_proportion, self.cap_percent, 0, self.complete, verbose=self.verbose, naive=False, intron_filter=self.intron_filter)
                    total_reads = locus.weight
                    if locus.graph:
                        locus.assemble_transcripts(complete=self.complete)
                    else:
                        continue
                    
                    if self.verbose:
                        reads_used = 0
                        transcripts_written = 0
                    
                    for transcript in locus.transcripts:
                        if self.passes_all_checks(transcript):
                            self.output_transcripts(transcript, self.output_type)
                            if self.verbose:
                                reads_used += transcript.weight
                                transcripts_written += 1
                    
                    if self.verbose:
                        print('\t{} transcripts assembled from {} of {} reads ({}%)'.format(
                            transcripts_written, round(reads_used,1), round(total_reads,1), round(reads_used/total_reads*100,2)))
                    
                    if subchunk[0].left() >= STOP_AT:
                        sys.exit()
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend assemble |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                                       {}\n".format(self.input)
        options_string += "  Output file (-o):                                 {}\n".format(self.output)
        options_string += "  Source name (--source):                           {}\n".format(self.source)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Max allowed gap in coverage (--max_gap):          {}\n".format(self.max_gap)
        options_string += "  Cluster distance for ends (--end_cluster):        {}\n".format(self.end_cluster)
        options_string += "  Min spanning bases (--min_overhang):              {}\n".format(self.min_overhang)
        options_string += "  Infer 5' ends from cov change (--infer_starts):   {}\n".format(self.infer_starts)
        options_string += "  Infer 3' ends from cov change (--infer_ends):     {}\n".format(self.infer_ends)
        options_string += "  *** Filters ***\n"
        options_string += "  Min bp transcript length (--minlen):              {}\n".format(self.minlen)
        options_string += "  Min isoform coverage (--min_cov):                 {}\n".format(self.min_cov)
        options_string += "  Min isoform contribution (--min_proportion):      {}\n".format(self.min_proportion)
        options_string += "  Min retained intron proportion (--intron_filter): {}\n".format(self.intron_filter)
        options_string += "  Min number of 5' reads (--min_start):             {}\n".format(self.min_start)
        options_string += "  Min number of 3' reads (--min_end):               {}\n".format(self.min_end)
        options_string += "  Min percent 5' end reads w/ uuG (--cap_percent):  {}\n".format(self.cap_percent)
        options_string += "  Keep fragmented assemblies (--allow_incomplete):  {}\n".format(self.incomplete)
        return options_string
    
    def display_summary(self):
        summary = ''
        summary += 'Total elapsed time: {}'.format(round(self.end_time - self.start_time, 5))
        return summary
    
    def file_extension(self, filename):
        """Boolean if the file's extension is valid (BED, ELR)"""
        split_name = filename.split('.')
        if len(split_name) == 1:
            return None
        else:
            return split_name[-1].lower()
    
    def input_is_valid(self, filename):
        """Boolean if the file is a format that Assembler can parse."""
        if self.file_extension(filename) in ['bed','elr','bam','sam','gtf','gff3','gff']:
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
        if transcript.get_length() < self.minlen: return False
        if transcript.attributes['S.reads'] < self.min_start: return False
        if transcript.attributes['E.reads'] < self.min_end: return False
        # if not args.INCOMPLETE and True not in transcript.splice and transcript.attributes['S.capped']/transcript.attributes['S.reads'] < args.CAP_PERCENT: return False
        if transcript.strand == '.' and transcript.coverage < self.min_unstranded: return False
        return True
    
    def run(self):
        """Executes end labeling on all reads."""
        print(self.display_options())
        for locus in self.generator:
            self.process_entry(locus)
        
        self.output_file.close()
        self.end_time = time.time()
        print(self.display_summary())

if __name__ == '__main__':
    from argument_parsers import assemble_parser as parser
    args = vars(parser.parse_args())
    obj = Assembler(args)
    sys.exit(obj.run())

