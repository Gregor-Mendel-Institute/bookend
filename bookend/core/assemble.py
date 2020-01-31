import re
import sys
import argparse
import time
import pysam
from cython_utils import _fasta_utils as fu
from cython_utils import _rnaseq_utils as ru
from cython_utils import _assembly_utils as au
if sys.version_info >= (3,0):
    izip = zip
else:
    from itertools import izip

###################
# INPUT ARGUMENTS #
###################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Assembles an end-labeled read BED12/ELR file. Produces an output assembly (BED12/ELR/GTF) and a table of summary statistics (bookend_stats.tsv)"
    )
    parser.add_argument(
        '-o','--output', dest='OUT',
        help="Destination file for assembly. File extension (bed, elr, gtf) determines output type.",
        type=str, default='bookend_assembly.gtf'
    )
    parser.add_argument(
        '--stats', dest='STATS',
        help="Destination file for summary statistics TSV file.",
        type=str, default='bookend_stats.tsv'
    )
    parser.add_argument(
        '--source', dest='SOURCE',
        help="String to write in the source column of output GTF/ELR",
        type=str, default='bookend'
    )
    parser.add_argument(
        '--genome', dest='GENOME',
        help="Genome FASTA file. Used for end label filtering only if input is BAM/SAM.",
        type=str, default=None
    )
    parser.add_argument(
        '--merge', dest='MERGE',
        help="Run assembly in merge mode.",
        default=False, action='store_true'
    )
    parser.add_argument(
        '--allow_incomplete', dest='INCOMPLETE',
        help="Filter out all assembled transcripts that are not end-to-end complete.",
        default=False, action='store_true'
    )
    parser.add_argument(
        '--max_gap', dest='MAX_GAP',
        help="Largest gap size to tolerate (number of nucleotides).",
        type=int, default=50
    )
    parser.add_argument(
        '--min_overhang', dest='MIN_OVERHANG',
        help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).",
        type=int, default=3
    )
    parser.add_argument(
        '--end_cluster', dest='END_CLUSTER',
        help="Largest distance between ends to consider in one cluster (number of nucleotides).",
        type=int, default=100
    )
    parser.add_argument(
        '--min_cov', dest='MIN_COV',
        help="Minimum coverage filter to remove low-evidence transcript models.",
        type=float, default=1.5
    )
    parser.add_argument(
        '--min_unstranded_cov', dest='MIN_UNSTRANDED',
        help="(Only used if --allow_incomplete) Set a more stringent threshold for keeping nonstranded frags.",
        type=float, default=20
    )
    parser.add_argument(
        '--min_start', dest='MIN_S',
        help="Minimum number of start reads to accept as a start site.",
        type=float, default=1
    )
    parser.add_argument(
        '--min_end', dest='MIN_E',
        help="Minimum number of end reads to accept as an end site",
        type=float, default=1
    )
    parser.add_argument(
        '--min_len', dest='MIN_LEN',
        help="Minimum output transcript length (nucleotides).",
        type=int, default=50
    )
    parser.add_argument(
        '--min_proportion', dest='MIN_PROPORTION',
        help="[float 0-1] Exclude ends, juctions, or transcripts that contribute < this proportion.",
        type=float, default=0.01
    )
    parser.add_argument(
        '--intron_filter', dest='INTRON_FILTER',
        help="[float 0-1] Retained introns must exceed this proportion the be considered.",
        type=float, default=0.15
    )
    parser.add_argument(
        '--cap_percent', dest='CAP_PERCENT',
        help="[float 0-1] Exclude 5' end features with < this percent of cap-labeled reads (reads that contain upstream untemplated G).",
        type=float, default=0.05
    )
    parser.add_argument(
        '--verbose', dest='VERBOSE',
        help="Display a verbose summary of each assembly in stdout.",
        default=False, action='store_true'
    )
    parser.add_argument(
        dest='INPUT',
        help="Input BED/ELR filepath. MUST be a single coordinate-sorted file.",
        type=str
    )
    args = parser.parse_args()
    return args

def file_extension(filename):
    """Boolean if the file's extension is valid (BED, ELR)"""
    split_name = filename.split('.')
    if len(split_name) == 1:
        return None
    else:
        return split_name[-1].lower()

def input_is_valid(filename):
    if file_extension(filename) in ['bed','elr','bam','sam']:
        return True
    else:
        return False

def passes_all_checks(transcript, args):
    """Determines if a transcript model passes the collection
    of filters set by the commandline arguments.
    'transcript' is an RNAseqMapping object,
    see _rnaseq_utils.pyx for details.
    """
    if transcript.coverage < args.MIN_COV: return False
    if not args.INCOMPLETE and not transcript.complete: return False
    if transcript.get_length() < args.MIN_LEN: return False
    if transcript.attributes['S.reads'] < args.MIN_S: return False
    if transcript.attributes['E.reads'] < args.MIN_E: return False
    # if not args.INCOMPLETE and True not in transcript.splice and transcript.attributes['S.capped']/transcript.attributes['S.reads'] < args.CAP_PERCENT: return False
    if transcript.strand == '.' and transcript.coverage < args.MIN_UNSTRANDED: return False
    return True

def output_transcripts(transcript, output_type, output):
    """Writes the RNAseqMapping object 'transcript' to an output stream,
    formatted as output_type."""
    if output_type == 'elr':
        output_line = transcript.write_as_elr()
    elif output_type == 'bed':
        output_line = transcript.write_as_bed(dataset.chrom_array, ['bookend'], score_column='coverage')
    elif output_type == 'gtf':
        output_line = transcript.write_as_gtf(dataset.chrom_array, 'bookend')
    
    output.write(output_line+'\n')

def output_statistics(transcript, stats, output):
    """Writes a TSV of the transcript's attributes to an output file connection."""
    line = '\t'.join([str(transcript.attributes.get(attr, 'NA')) for attr in stats])+'\n'
    output.write(line)

################
# STATS CONFIG #
################

stats_config = [
    'transcript_id',
    'gene_id',
    'cov',
    'length',
    'reads',
    'S.reads',
    'S.capped',
    'S.left',
    'S.right',
    'E.reads',
    'E.left',
    'E.right'
]

###################################

if __name__ == '__main__':
    args = parse_args()
    if args.INCOMPLETE:
        # Some settings are incompatible with writing incomplete transcripts.
        # Adjust args to allow incomplete.
        args.MIN_S = 0
        args.MIN_E = 0
    
    if input_is_valid(args.INPUT):
        file_type = file_extension(args.INPUT)
        if file_type in ['bam','sam']:
            input_file = pysam.AlignmentFile(args.INPUT)
            dataset = ru.RNAseqDataset(genome_fasta=args.GENOME, chrom_array=input_file.header.references)
        else:
            dataset = ru.RNAseqDataset()
            input_file = open(args.INPUT)
    else:
        print("\nERROR: input file must be a valid format (BED, ELR, BAM, SAM).")
        parser.print_help()
        sys.exit(1)
    
    stats_output = open(args.STATS, 'w')
    stats_output.write('\t'.join(stats_config)+'\n')
    STOP_AT=float('inf')
    complete = not args.INCOMPLETE
    # STOP_AT=62485
    # STOP_AT=1911044
    output_type = file_extension(args.OUT)
    if output_type is None:
        output_type = 'gtf'
    
    start_time = time.time()
    chunk_counter = 0
    chunks = ru.read_generator(input_file, dataset, file_type, args.MAX_GAP)
    if args.VERBOSE:print("Initialized generator ({})".format(time.time()-start_time))
    output = open(args.OUT,'w')
    for chunk in chunks:
        if len(chunk) > 0:
            chrom = chunk[0].chrom
            subchunks = ru.split_chunk(chunk, args.MIN_PROPORTION, args.MAX_GAP)
            for subchunk in subchunks:
                chunk_counter += 1
                if len(subchunk) > 0:
                    if args.VERBOSE:print('\t{}\n[{}:{}-{}] Processing chunk.'.format(time.time()-start_time, dataset.chrom_array[chrom], subchunk[0].left(),subchunk[-1].right()))
                    start_time = time.time()
                    locus = au.Locus(chrom, chunk_counter, subchunk, args.MAX_GAP, args.END_CLUSTER, args.MIN_OVERHANG, True, args.MIN_PROPORTION, args.CAP_PERCENT, 0, complete, verbose=args.VERBOSE, naive=False, intron_filter=args.INTRON_FILTER)
                    total_reads = locus.weight
                    if locus.graph:
                        locus.assemble_transcripts(complete=complete)
                    else:
                        continue
                    
                    if args.VERBOSE:
                        reads_used = 0
                        transcripts_written = 0
                    
                    for transcript in locus.transcripts:
                        if passes_all_checks(transcript, args):
                            output_transcripts(transcript, output_type, output)
                            output_statistics(transcript, stats_config, stats_output)
                            if args.VERBOSE:
                                reads_used += transcript.weight
                                transcripts_written += 1
                    
                    if args.VERBOSE:
                        print('\t{} transcripts assembled from {} of {} reads ({}%)'.format(
                            transcripts_written, round(reads_used,1), round(total_reads,1), round(reads_used/total_reads*100,2)))
                    
                    if subchunk[0].left() >= STOP_AT:
                        sys.exit()
    
    output.close()

