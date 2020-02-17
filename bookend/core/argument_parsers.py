import argparse
from argparse import RawTextHelpFormatter
from bookend.core.assemble import Assembler
from bookend.core.fastq_end_label import EndLabeler
from bookend.core.bam_to_elr import BAMtoELRconverter
from bookend.core.fasta_index import Indexer
from bookend.core.bed_to_elr import BEDtoELRconverter
from bookend.core.elr_to_bed import ELRtoBEDconverter
# from core.elr_combine import ELRcombiner
# from core.sam_sj_out import SAMtoSJconverter
# from core.sj_merge import SJmerger
# from core.sj_to_bed import SJtoBEDconverter
# from bookend.core.gtf_merge import AnnotationMerger

#TODO: Flesh out Helper class
class Helper:
    def __init__(self, args):
        self.no_args_message = """
/| bookend |\\
¯¯¯¯¯¯¯¯¯¯¯¯¯
End-guided transcriptome assembly.
usage: bookend [subcommand] [options] [input file(s)]
Subcommands (use -h/--help for more info):

    label    (Label 5' and 3' ends in a FASTQ file)
    assemble (Assemble transcripts from aligned end-labeled reads)
    merge    (Merge assembled GTFs with a reference annotation)

    --end-labeled read (ELR) operations--
    make-elr
    sort-elr
    collapse-elr
    combine-elr

    --file conversion--
    bed-to-elr
    elr-to-bed
    gff3-to-bed
    gtf-to-bed
    sam-to-sj
    sj-to-bed
    merge-sj

    --setup/indexing--
    index-fasta
    softbridge-fasta

"""
    
    def run(self):
        print(self.no_args_message)


desc = 'Functions for end-guided assembly of RNA-seq data.'
main_parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
main_parser.set_defaults(object=Helper)
subparsers = main_parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Supported subcommands:')

ELRdesc = """
Converts file to a 'end labeled read' (ELR) format.
ELR contains a two-component header: #C (chromosome) and #S (source).
Each line is a read or read stack with six columns:
    chromosome  position  strand  ELCIGAR  sample  weight

    chromosome: Chromosome number, as indexed by the #C header
    position: 0-indexed position in chromosome
    strand: Inferred RNA strand; +, -, or . (unknown)
    ELCIGAR: String describing the mapped read (full description below)
    sample: Sample number, as indexed by the #S header
    weight: Read count (allows partial counts for multimappers)
"""

ELRepilog = """
ELCIGAR strings are Character|Number strings and one trailing character
([CN]xC), where C is a label and N is a numeric length on the genome.
Each C labels an end or gap in the alignment as one of the following:
    S: start
    s: start (low confidence)
    C: start (capped)
    E: end (polyA tail)
    e: end (low confidence)
    D: splice junction donor
    A: splice junction acceptor
    .: unspecified gap or end

For example, a 50bp paired-end read of a 185-nt fragment would be
    .50.85.50.
A full 3-exon transcript could be described as:
    S256D800A128D800A512E
where the 3 exons are 256, 128, and 512 nucleotides,
and the 2 introns are both 800 nucleotides.
"""

### assemble.py ###
assemble_parser = subparsers.add_parser('assemble',help="Assembles an end-labeled read (ELR) file. Produces an output assembly (BED12/ELR/GTF) and a table of summary statistics (bookend_stats.tsv)")
assemble_parser.add_argument('-o','--output', dest='OUT', type=str, default='bookend_assembly.gtf', help="Destination file for assembly. File extension (bed, elr, gtf) determines output type.")
assemble_parser.add_argument('--source', dest='SOURCE', type=str, default='bookend', help="String to write in the source column of output GTF/ELR")
assemble_parser.add_argument('--genome', dest='GENOME', type=str, default=None, help="Genome FASTA file. Used for end label filtering only if input is BAM/SAM.")
assemble_parser.add_argument('--max_gap', dest='MAX_GAP', type=int, default=50, help="Largest gap size to tolerate (number of nucleotides).")
assemble_parser.add_argument('--min_overhang', dest='MIN_OVERHANG', type=int, default=3, help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).")
assemble_parser.add_argument('--allow_incomplete', dest='INCOMPLETE', default=False, action='store_true', help="Filter out all assembled transcripts that are not end-to-end complete.")
assemble_parser.add_argument('--infer_starts', dest='INFER_STARTS', default=False, action='store_true', help="If S tags are missing, calculate the most likely start site based on coverage changes.")
assemble_parser.add_argument('--infer_ends', dest='INFER_ENDS', default=False, action='store_true', help="If E tags are missing, calculate the most likely end site based on coverage changes.")
assemble_parser.add_argument('--end_cluster', dest='END_CLUSTER', type=int, default=100, help="Largest distance between ends to consider in one cluster (number of nucleotides).")
assemble_parser.add_argument('--min_cov', dest='MIN_COV', type=float, default=1.5, help="Minimum coverage filter to remove low-evidence transcript models.")
assemble_parser.add_argument('--min_unstranded_cov', dest='MIN_UNSTRANDED', type=float, default=20, help="(Only used if --allow_incomplete) Set a more stringent threshold for keeping nonstranded frags.")
assemble_parser.add_argument('--min_start', dest='MIN_S', type=float, default=1, help="Minimum number of start reads to accept as a start site.")
assemble_parser.add_argument('--min_end', dest='MIN_E', type=float, default=1, help="Minimum number of end reads to accept as an end site")
assemble_parser.add_argument('--minlen', dest='MINLEN', type=int, default=50, help="Minimum output transcript length (nucleotides).")
assemble_parser.add_argument('--min_proportion', dest='MIN_PROPORTION', type=float, default=0.01, help="[float 0-1] Exclude ends, juctions, or transcripts that contribute < this proportion. (Used as a signal threshold)")
assemble_parser.add_argument('--intron_filter', dest='INTRON_FILTER', type=float, default=0.15, help="[float 0-1] Retained introns must exceed this proportion the be considered.")
assemble_parser.add_argument('--cap_percent', dest='CAP_PERCENT', type=float, default=0.05, help="[float 0-1] Exclude 5' end features with < this percent of cap-labeled reads (reads that contain upstream untemplated G).")
assemble_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each assembly in stdout.")
assemble_parser.add_argument(dest='INPUT', type=str, help="Input BED/ELR filepath. MUST be a single coordinate-sorted file.")
assemble_parser.set_defaults(object=Assembler)

### bam_to_elr.py ###
bam_to_elr_parser = subparsers.add_parser('make-elr',help="Converts a BAM or SAM file to an End-Labeled Read (ELR) or BED12 file.", description=ELRdesc, epilog=ELRepilog, formatter_class=RawTextHelpFormatter)
bam_to_elr_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, default=None, help="Filepath to write end-labeled file.")
bam_to_elr_parser.add_argument("--source", dest='SOURCE', default=None, type=str, help="Name the source of BAM/SAM reads.")
bam_to_elr_parser.add_argument("--genome", dest='GENOME', default=None, type=str, help="Genome FASTA file")
bam_to_elr_parser.add_argument("--stranded", dest='STRANDED', default=False, action="store_true", help="The reads are strand-specific")
bam_to_elr_parser.add_argument("-s", dest='START', default=False, action='store_true', help="Read 5' ends are transcript start sites.")
bam_to_elr_parser.add_argument("-c", dest='CAPPED', default=False, action='store_true', help="5' end data is capped.")
bam_to_elr_parser.add_argument("-e", dest='END', default=False, action='store_true', help="Read 3' ends are transcript end sites.")
bam_to_elr_parser.add_argument("--no_ends", dest='NO_ENDS', default=False, action='store_true', help="(overrides other options) Ignore all 5' and 3' end labels.")
bam_to_elr_parser.add_argument("--bed", dest='BED_OUT', default=False, action='store_true', help="Output a 15-column end labeled BED file.")
bam_to_elr_parser.add_argument("--secondary", dest='SECONDARY', default=False, action='store_true', help="Keep secondary alignments.")
bam_to_elr_parser.add_argument("--header", dest='HEADER', type=str, default=None, help="Filepath to write ELR header.")
bam_to_elr_parser.add_argument("--start_seq", dest='START_SEQ', default='ACGGG', type=str, help="Sequence of the oligo that marks a 5' read (sense)")
bam_to_elr_parser.add_argument("--end_seq", dest='END_SEQ', default='RRRRRRRRRRRRRRR', type=str, help="Sequence of the oligo that marks a 3' read (sense)")
bam_to_elr_parser.add_argument("--record_artifacts", dest='RECORD_ARTIFACTS', default=False, action='store_true', help="Reports artifact-masked S/E labels as lowercase s/e.")
bam_to_elr_parser.add_argument("--split", dest='SPLIT', default=False, action='store_true', help="Separate reads into different files by their multimapping number.")
bam_to_elr_parser.add_argument("--mismatch_rate", dest='MM_RATE', default=0.20, type=float, help="Mismatch tolerance of S label matches")
bam_to_elr_parser.add_argument("--sj_shift", dest='SJ_SHIFT', default=2, type=int, help="Shift up to this many bases to find a canonical splice junction")
bam_to_elr_parser.add_argument("--minlen", dest='MINLEN', default=20, type=int, help="Minimum number of nucleotides for an aligned read.")
bam_to_elr_parser.add_argument("INPUT", type=str, default=None, help="Input BAM/SAM file")
bam_to_elr_parser.set_defaults(object=BAMtoELRconverter)

### bed_to_elr.py ###
bed_to_elr_parser = subparsers.add_parser('bed-to-elr',help="Converts a BED file to an End-Labeled Read (ELR) file.", description=ELRdesc, epilog=ELRepilog, formatter_class=RawTextHelpFormatter)
bed_to_elr_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, required=True, help="Filepath to write ELR file.")
bed_to_elr_parser.add_argument("--header", dest='HEADER', type=str, default=None, help="Filepath to write ELR header.")
bed_to_elr_parser.add_argument("--source", dest='SOURCE', help="Source of BED lines.", default=None, type=str)
bed_to_elr_parser.add_argument("-j", dest='JUNCTIONS', default=False, action='store_true', help="Gaps in the reads are all splice junctions.")
bed_to_elr_parser.add_argument("-s", dest='START', default=False, action='store_true', help="All read 5' ends are transcript start sites.")
bed_to_elr_parser.add_argument("-c", dest='CAPPED', default=False, action='store_true', help="All 5' ends are capped.")
bed_to_elr_parser.add_argument("-e", dest='END', default=False, action='store_true', help="All read 3' ends are transcript end sites.")
bed_to_elr_parser.add_argument("INPUT", type=str, help='Input BED file')
bed_to_elr_parser.set_defaults(object=BEDtoELRconverter)

### elr_combine.py ###

# parser.add_argument(dest='INPUT', help="Input sorted ELR file.", type=str, nargs='+')

# args = parser.parse_args()

### elr_to_bed.py ###
elr_to_bed_parser = subparsers.add_parser('elr-to-bed',help="Converts an End-Labeled Read (ELR) file to BED12.", formatter_class=RawTextHelpFormatter)
elr_to_bed_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, default=None, required=True, help="Filepath to write BED file.")
elr_to_bed_parser.add_argument("--header", dest='HEADER', type=str, default=None, help="Filepath to write ELR header.")
elr_to_bed_parser.add_argument("INPUT", help='Input ELR file')
elr_to_bed_parser.set_defaults(object=ELRtoBEDconverter)

### fastq_end_label.py ###    
end_label_parser = subparsers.add_parser('label',help="Trims and labels RNA 5' and 3' ends in a FASTQ file")
end_label_parser.add_argument('-S', '--start', dest='START', type=str, default='AAGCAGTGGTATCAACGCAGAGTACGGG', help="Template switching primer, marks the RNA 5' terminus.")
end_label_parser.add_argument( '-E', '--end', dest='END', type=str, default='AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTT+', help="Reverse transcription primer, marks (reverse complement of) the RNA 3' terminus.")
end_label_parser.add_argument('--strand', dest='STRAND', type=str, default='unstranded', choices=['forward','reverse','unstranded'], help="Orientation of mate1 with respect to the RNA strand")
end_label_parser.add_argument('--out1', dest='OUT1', type=str, default='end_label.1.fastq', help="Destination file for trimmed and labeled mate 1 FASTQ file.")
end_label_parser.add_argument('--out2', dest='OUT2', type=str, default='end_label.2.fastq', help="Destination file for trimmed and labeled mate 2 FASTQ file.")
end_label_parser.add_argument('--single_out', dest='SINGLE_OUT', type=str, default='end_label.single.fastq', help="Destination file for reads that are only informative single-end.")
end_label_parser.add_argument('--min_start', dest='MIN_START', type=int, default=7, help="Minimum number of nucleotides that must match to the start adapter sequence.")
end_label_parser.add_argument('--min_end', dest='MIN_END', type=int, default=9, help="Minimum number of nucleotides that must match to the end adapter sequence.")
end_label_parser.add_argument('--suppress_untrimmed', dest='SUPPRESS_UNTRIMMED', default=False, action='store_true', help="If no trimming occurred, do not write the read to output.")
end_label_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display each trimming result on stdout.")
end_label_parser.add_argument('--minlen', dest='MINLEN', type=int, default=16, help="Minimum sequence length to keep.")
end_label_parser.add_argument('--mismatch_rate', dest='MM_RATE', type=float, default=.06, help="Highest allow proportion of mismatches.")
end_label_parser.add_argument('--qualmask', dest='QUALMASK', type=float, default=16, help="Ignores any basecalls with phred score < this, treats base as N.")
end_label_parser.add_argument('--minqual', dest='MINQUAL', type=float, default=25, help="Suppresses any trimmed sequences with lower than this mean phred quality score.")
end_label_parser.add_argument(dest='FASTQ', type=str, nargs='+', help="Input FASTQ file(s). 1 for single-end, 2 for paired-end")
end_label_parser.set_defaults(object=EndLabeler)

### fasta_index.py ###
fasta_index_parser = subparsers.add_parser('index-fasta',help="Generates index files for the reference genome FASTA file.")
fasta_index_parser.add_argument('FASTA', type=str,help="Path to FASTA file.")
fasta_index_parser.add_argument('-S','--split', dest='SPLIT_ON', type=str, help="Character string to split header name from description.",default=' ')
fasta_index_parser.add_argument("--bridge_min", dest='MINLEN',help="Minimum tolerated length for a bridge over a softmasked region.", default=100, type=int)
fasta_index_parser.add_argument("--bridge_max", dest='MAXLEN',help="Maximum tolerated length for a bridge over a softmasked region.", default=500, type=int)
fasta_index_parser.set_defaults(object=Indexer)

### gtf_merge.py ###
# merge_parser = subparsers.add_parser('merge',help="Merges multiple assembly/annotation GTF/GFF3 files to one consensus")
# merge_parser.add_argument('-r', dest='reference_GFF', help="[GFF3/GTF] Path to reference annotation(s)", nargs='*')
# merge_parser.add_argument('-i', dest='input_GFF', help="[GFF3/GTF] Path to input annotation(s)", nargs='+')
# merge_parser.add_argument('--genome', dest='genome_fasta', help="Path to a FASTA file of the genome.", default=None)
# merge_parser.add_argument("--tracking", dest='tracking_file', help="gffcompare .tracking file.", default=None)
# merge_parser.add_argument("--stats", dest='stats_file', help="bookend .stats file.", default=None)
# merge_parser.add_argument("--buffer", dest='buffer', help="Number of allowed nonoverlapping terminal nucleotides.", default=100, type=int)
# merge_parser.add_argument("--fasta_out", dest='fasta_out', help="Output FASTA file of transcript sequences.", default=None, type=str)
# merge_parser.add_argument("--protein_out", dest='protein_out', help="Output FASTA file of amino acid sequences.", default=None, type=str)
# merge_parser.set_defaults(object=AnnotationMerger)

### sam_sj_out.py ###
# parser.add_argument(
#     "-F", "--fasta", dest='FASTA',
#     help="Genome FASTA file",
#     default=None, type=str, required=True
# )
# parser.add_argument(
#     "--format", dest='FORMAT',
#     help="Output file format",
#     default='star', type=str, choices=['bed','star']
# )
# parser.add_argument(
#     "--filter", dest='FILTER',
#     help="Remove noncanonical splice junctions from the output",
#     default=False, action='store_true'
# )
# parser.add_argument(
#     "FILENAME", nargs='?'
# )

### sj_merge.py ###

### sj_to_bed.py ###


