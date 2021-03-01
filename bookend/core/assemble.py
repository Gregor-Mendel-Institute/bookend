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


# assemble_parser = subparsers.add_parser('assemble',help="Assembles an end-labeled read (ELR) file. Produces an output assembly (BED12/ELR/GTF) and a table of summary statistics (bookend_stats.tsv)")
# assemble_parser.add_argument('-o','--output', dest='OUT', type=str, default='bookend_assembly.gtf', help="Destination file for assembly. File extension (bed, elr, gtf) determines output type.")
# assemble_parser.add_argument('--source', dest='SOURCE', type=str, default='bookend', help="String to write in the source column of output GTF/ELR")
# assemble_parser.add_argument('--genome', dest='GENOME', type=str, default=None, help="Genome FASTA file. Used for end label filtering only if input is BAM/SAM.")
# assemble_parser.add_argument('--max_gap', dest='MAX_GAP', type=int, default=50, help="Largest gap size to tolerate (number of nucleotides).")
# assemble_parser.add_argument('--min_overhang', dest='MIN_OVERHANG', type=int, default=3, help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).")
# assemble_parser.add_argument('--allow_incomplete', dest='INCOMPLETE', default=False, action='store_true', help="Keep assembled transcripts even if they are not end-to-end complete.")
# assemble_parser.add_argument('--infer_starts', dest='INFER_STARTS', default=False, action='store_true', help="If S tags are missing, calculate the most likely start site based on coverage changes.")
# assemble_parser.add_argument('--infer_ends', dest='INFER_ENDS', default=False, action='store_true', help="If E tags are missing, calculate the most likely end site based on coverage changes.")
# assemble_parser.add_argument('--end_cluster', dest='END_CLUSTER', type=int, default=100, help="Largest distance between ends to consider in one cluster (number of nucleotides).")
# assemble_parser.add_argument('--min_cov', dest='MIN_COV', type=float, default=1.5, help="Minimum coverage filter to remove low-evidence transcript models.")
# assemble_parser.add_argument('--min_unstranded_cov', dest='MIN_UNSTRANDED', type=float, default=20, help="(Only used if --allow_incomplete) Set a more stringent threshold for keeping nonstranded frags.")
# assemble_parser.add_argument('--min_start', dest='MIN_S', type=float, default=1, help="Minimum number of start reads to accept as a start site.")
# assemble_parser.add_argument('--min_end', dest='MIN_E', type=float, default=1, help="Minimum number of end reads to accept as an end site")
# assemble_parser.add_argument('--minlen', dest='MINLEN', type=int, default=50, help="Minimum output transcript length (nucleotides).")
# assemble_parser.add_argument('--min_proportion', dest='MIN_PROPORTION', type=float, default=0.01, help="[float 0-1] Exclude ends, juctions, or transcripts that contribute < this proportion. (Used as a signal threshold)")
# assemble_parser.add_argument('--intron_filter', dest='INTRON_FILTER', type=float, default=0.15, help="[float 0-1] Retained introns must exceed this proportion the be considered.")
# assemble_parser.add_argument('--cap_bonus', dest='CAP_BONUS', type=float, default=5, help="[float] Signal multiplier for 5' reads with an inferred cap structure (uuG).")
# assemble_parser.add_argument("--ignore_labels", dest='IGNORE_LABELS', default=False, action='store_true', help="(overrides other options) Ignore all 5' and 3' end labels.")
# assemble_parser.add_argument("--ignore_source", dest='IGNORE_SOURCE', default=False, action='store_true', help="Do not separate read weights by source.")
# assemble_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each assembly in stdout.")
# assemble_parser.add_argument(dest='INPUT', type=str, help="Input BED/ELR filepath. MUST be a single coordinate-sorted file.")
# assemble_parser.set_defaults(object='Assembler')

class Assembler:
    def __init__(self, args):
        """Parses input arguments for assembly"""
        self.start_time = time.time()
        self.output = args['OUT']
        self.source = args['SOURCE']
        self.genome = args['GENOME']
        self.incomplete = args['INCOMPLETE']
        self.max_gap = args['MAX_GAP']
        self.end_cluster = args['END_CLUSTER']
        self.min_overhang = args['MIN_OVERHANG']
        self.min_cov = args['MIN_COV']
        self.min_unstranded_cov = args['MIN_UNSTRANDED']
        self.min_start = args['MIN_S']
        self.min_end = args['MIN_E']
        self.intron_filter = args['INTRON_FILTER']
        self.min_proportion = args['MIN_PROPORTION']
        self.cap_bonus = args['CAP_BONUS']
        self.minlen = args['MINLEN']
        self.verbose = args['VERBOSE']
        self.input = args['INPUT']
        self.ignore_labels = args['IGNORE_LABELS']
        self.ignore_source = args['IGNORE_SOURCE']
        
        self.naive = not self.ignore_source
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
        
        self.generator = ru.read_generator(self.input_file, self.dataset, self.file_type, self.max_gap, self.min_proportion)
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
            self.chunk_counter += 1
            if self.verbose:
                print('\n[{}:{}-{}] Processing chunk.'.format(self.dataset.chrom_array[chrom], chunk[0].left(),chunk[-1].right()))
            
            locus = au.Locus(chrom, self.chunk_counter, chunk, self.max_gap, self.end_cluster, self.min_overhang, True, self.min_proportion, self.cap_bonus, self.complete, verbose=self.verbose, naive=self.naive, intron_filter=self.intron_filter, ignore_ends=self.ignore_labels)
            total_reads = locus.weight
            if locus.graph:
                locus.assemble_transcripts(complete=self.complete)
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
                
                if chunk[0].left() >= STOP_AT:
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
        options_string += "  Ignore read sources (--ignore_source):            {}\n".format(self.ignore_source)
        options_string += "  Ignore end labels (--ignore_labels):              {}\n".format(self.ignore_labels)
        options_string += "  *** Filters ***\n"
        options_string += "  Min bp transcript length (--minlen):              {}\n".format(self.minlen)
        options_string += "  Min isoform coverage (--min_cov):                 {}\n".format(self.min_cov)
        options_string += "  Min isoform contribution (--min_proportion):      {}\n".format(self.min_proportion)
        options_string += "  Min retained intron proportion (--intron_filter): {}\n".format(self.intron_filter)
        options_string += "  Min number of 5' reads (--min_start):             {}\n".format(self.min_start)
        options_string += "  Min number of 3' reads (--min_end):               {}\n".format(self.min_end)
        options_string += "  Min percent 5' end reads w/ uuG (--cap_bonus):  {}\n".format(self.cap_bonus)
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

