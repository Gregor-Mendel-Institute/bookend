#!/usr/bin/python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser, RawTextHelpFormatter, ArgumentDefaultsHelpFormatter

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
    elr      (Convert a BAM/SAM file to the end-labeled read format)
    assemble (Assemble transcripts from aligned end-labeled reads)
    condense (Partial assembly that leaves keeps all fragments; use for meta-assembly)
    classify (Compare an assembly to the transcripts of a reference annotation)
    merge    (Unify a set of assemblies with or without a reference annotation)
    quantify (Count the abundance of assembled transcripts in a given ELR/BAM file)
    bedgraph (Write a coverage Bedgraph file of end-labeled reads)
    fasta    (Write a transcript FASTA file from an annotation and genome)

    --end-labeled read (ELR) operations--
    elr-sort
    elr-subset
    elr-combine
    
    --file conversion--
    gtf-to-bed
    bed-to-elr
    elr-to-bed
    sam-to-sj
    sj-to-bed
    sj-merge

"""
    
    def run(self):
        print(self.no_args_message)


desc = 'Functions for end-guided assembly of RNA-seq data.'
main_parser = ArgumentParser(description=desc, formatter_class=ArgumentDefaultsHelpFormatter)
main_parser.set_defaults(object='Helper')
subparsers = main_parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Supported subcommands:')

ELRdesc = """
Converts file to a 'end labeled read' (ELR) format.
ELR contains a two-component header: #C (chromosome) and #S (source).
Each line is a read or read stack with seven columns:
    chromosome  position  length strand  ELCIGAR  sample  weight

    chromosome: Chromosome number, as indexed by the #C header
    position: 0-indexed position in chromosome
    length: Distance between the leftmost and righmost mapped position
    strand: Inferred RNA strand; +, -, or . (unknown)
    ELCIGAR: String describing the mapped read (full description below)
    sample: Sample number, as indexed by the #S header
    weight: Read count (allows partial counts for multimappers)
"""


### assemble.py ###
assemble_parser = subparsers.add_parser('assemble',help="Assembles an end-labeled read (ELR) file. Produces an output assembly (BED12/ELR/GTF) and a table of summary statistics (bookend_stats.tsv)", formatter_class=ArgumentDefaultsHelpFormatter)
assemble_parser.add_argument('-o','--output', dest='OUT', type=str, default='bookend_assembly.gtf', help="Destination file for assembly. File extension (bed, elr, gtf) determines output type.")
assemble_parser.add_argument("--source", dest='SOURCE', default='bookend', type=str, help="Name to add to the GTF source column.")
assemble_parser.add_argument('--cov_out', dest='COV_OUT', type=str, default=None, help="Destination for a TSV of coverage estimates for each transcript in each source.")
assemble_parser.add_argument('--max_gap', dest='MAX_GAP', type=int, default=50, help="Largest gap size to tolerate (nucleotides).")
assemble_parser.add_argument('--max_intron', dest='MAX_INTRON', type=int, default=100000, help="Ignore reads with introns longer than this.")
assemble_parser.add_argument('--end_cluster', dest='END_CLUSTER', type=int, default=100, help="Largest distance between end-labeled reads to consider the same cluster (nucleotides).")
assemble_parser.add_argument('--min_overhang', dest='MIN_OVERHANG', type=int, default=3, help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).")
assemble_parser.add_argument('--min_cov', dest='MIN_COV', type=float, default=2, help="Minimum coverage filter to remove low-evidence transcript models.")
assemble_parser.add_argument('--allow_incomplete', dest='INCOMPLETE', default=False, action='store_true', help="Keep assembled transcripts even if they are not end-to-end complete.")
assemble_parser.add_argument('--min_unstranded_cov', dest='MIN_UNSTRANDED', type=float, default=20, help="(Only used if --allow_incomplete) Set a more stringent threshold for keeping nonstranded frags.")
assemble_parser.add_argument('--min_start', dest='MIN_S', type=float, default=0, help="Minimum number of start reads to accept as a start site.")
assemble_parser.add_argument('--min_end', dest='MIN_E', type=float, default=0, help="Minimum number of end reads to accept as an end site")
assemble_parser.add_argument('--min_intron_len', dest='MIN_INTRON_LEN', type=int, default=50, help="Minimum length for an intron (nucleotides)")
assemble_parser.add_argument('--min_len', dest='MINLEN', type=int, default=100, help="Minimum output transcript length (nucleotides).")
assemble_parser.add_argument('--min_proportion', dest='MIN_PROPORTION', type=float, default=0.01, help="[float 0-1] Exclude ends, juctions, or transcripts that contribute < this proportion. (Used as a signal threshold)")
assemble_parser.add_argument('--intron_filter', dest='INTRON_FILTER', type=float, default=0.15, help="[float 0-1] Retained introns must exceed this proportion the be considered.")
assemble_parser.add_argument('--truncation_filter', dest='TRUNCATION_FILTER', type=float, default=0.5, help="[float 0-1] 5'-truncated transcripts must exceed this proportion of coverage of longer transcripts. 1 suppresses truncations, 0 allows all.")
assemble_parser.add_argument('--cap_bonus', dest='CAP_BONUS', type=float, default=5, help="[float] Signal multiplier for 5' reads with an inferred cap structure (uuG).")
assemble_parser.add_argument('--cap_filter', dest='CAP_FILTER', type=float, default=.02, help="[float] Require putative truncations to have >= this percent 'capped' reads (uuG).")
assemble_parser.add_argument("--use_sources", dest='USE_SOURCES', default=False, action='store_true', help="Attempt to maximize coherence of reads from the same source.")
assemble_parser.add_argument("--ignore_labels", dest='IGNORE_LABELS', default=False, action='store_true', help="(overrides other options) Ignore all 5' and 3' end labels.")
assemble_parser.add_argument("--require_cap", dest='REQUIRE_CAP', default=False, action='store_true', help="No start site is allowed to have less than cap_filter of uuG reads.")
assemble_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each assembly in stdout.")
assemble_parser.add_argument(dest='INPUT', type=str, nargs='+', help="Input ELR filepath(s). MUST be coordinate-sorted.")
assemble_parser.set_defaults(object='Assembler')

### bedgraph.py ###
bedgraph_parser = subparsers.add_parser('bedgraph',help="Produces a Bedgraph file from an end-labeled read (ELR) file.", formatter_class=ArgumentDefaultsHelpFormatter)
bedgraph_parser.add_argument('-o','--output', dest='OUT', type=str, default='bookend.bedgraph', help="Bedgraph destination file.")
bedgraph_parser.add_argument('-t','--type', dest='TYPE', choices=['','COV','5P','3P','S', 'E', 'C'], type=str, default='COV', help="Coverage type: filters Bedgraph output by read label. Default all")
bedgraph_parser.add_argument('-s','--strand', dest='STRAND', type=str, choices=['.','+','-'], default='.', help="Strand (-|+|.), default . (unstranded)")
bedgraph_parser.add_argument('--infer', dest='INFER_STRAND', default=False, action="store_true", help="Infer the strand of unstranded reads based on the strand ratio of the locus.")
bedgraph_parser.add_argument('--offset', dest='OFFSET', type=int, default=0, help="Bases downstream (+) or upstream (-) to offset 5P or 3P position.")
# bedgraph_parser.add_argument('--scale', dest='SCALE', default=False, action='store_true', help="Perform per-million scaling on output values.")
bedgraph_parser.add_argument(dest='INPUT', type=str, nargs='+', help="Input ELR filepath(s). MUST be coordinate-sorted.")
bedgraph_parser.set_defaults(object='Bedgrapher')

### bam_to_elr.py ###
bam_to_elr_parser = subparsers.add_parser('elr',help="Converts a BAM or SAM file to an End-Labeled Read (ELR) or BED12 file.", description=ELRdesc, formatter_class=ArgumentDefaultsHelpFormatter)
bam_to_elr_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, default=None, help="Filepath to write end-labeled file.")
bam_to_elr_parser.add_argument("--data_type", dest='DATA_TYPE', default=None, type=str, help="Choose presets for certain end-labeled libraries (smartseq, pacbio, ont, direct_rna)")
bam_to_elr_parser.add_argument("--source", dest='SOURCE', default=None, type=str, help="Name the source of BAM/SAM reads.")
bam_to_elr_parser.add_argument("--genome", dest='GENOME', default=None, type=str, help="Genome FASTA file")
bam_to_elr_parser.add_argument('--reference', dest='REFERENCE', help="[GFF3/GTF] Path to reference annotation", type=str, default=None)
bam_to_elr_parser.add_argument("--splice", dest='SPLICE', default=None, type=str, help="Reference splice junction file (STAR SJ.out.tab, intron BED6)")
bam_to_elr_parser.add_argument('--max_intron', dest='MAX_INTRON', type=int, default=100000, help="Ignore reads with introns longer than this.")
bam_to_elr_parser.add_argument('--max_indel', dest='MAX_INDEL', type=int, default=10, help="Ignore reads with indels longer than this.")
bam_to_elr_parser.add_argument("--chrom_names", dest='CHROM_NAMES', default=None, type=str, help="2-column TSV of chromosome ID -> chromosome name.")
bam_to_elr_parser.add_argument("--start_seq", dest='START_SEQ', default='ACATGGG', type=str, help="Sequence of the oligo that marks a 5' read (sense)")
bam_to_elr_parser.add_argument("--end_seq", dest='END_SEQ', default='RRRRRRRRRRRRRRRRRRRRRRRRRRRRRR', type=str, help="Sequence of the oligo that marks a 3' read (sense)")
bam_to_elr_parser.add_argument("--mismatch_rate", dest='MM_RATE', default=0.20, type=float, help="Mismatch tolerance of S label matches")
bam_to_elr_parser.add_argument("--stranded", dest='STRANDED', default=False, action="store_true", help="The reads are strand-specific")
bam_to_elr_parser.add_argument("--reverse", dest='REVERSE', default=False, action="store_true", help="Reads are oriented mate1-reverse")
bam_to_elr_parser.add_argument("-s", dest='START', default=False, action='store_true', help="Read 5' ends are transcript start sites.")
bam_to_elr_parser.add_argument("-c", dest='CAPPED', default=False, action='store_true', help="5' end data is capped.")
bam_to_elr_parser.add_argument("-e", dest='END', default=False, action='store_true', help="Read 3' ends are transcript end sites.")
bam_to_elr_parser.add_argument("--no_ends", dest='NO_ENDS', default=False, action='store_true', help="(overrides other options) Ignore all 5' and 3' end labels.")
bam_to_elr_parser.add_argument("--secondary", dest='SECONDARY', default=False, action='store_true', help="Keep secondary alignments.")
bam_to_elr_parser.add_argument("--header", dest='HEADER', type=str, default=None, help="Filepath to write ELR header.")
bam_to_elr_parser.add_argument("--record_artifacts", dest='RECORD_ARTIFACTS', default=False, action='store_true', help="Reports artifact-masked S/E labels as >/].")
bam_to_elr_parser.add_argument("--split", dest='SPLIT', default=False, action='store_true', help="Separate reads into different files by their multimapping number.")
bam_to_elr_parser.add_argument("--untrimmed", dest='UNTRIMMED', default=False, action='store_true', help="(overrides -s -c -e) End labels were not trimmed from input reads prior to alignment.")
bam_to_elr_parser.add_argument("--allow_noncanonical", dest='ALLOW_NONCANONICAL', default=False, action='store_true', help="Do not require canonical splice junction motifs (GT/AG, GC/AG, AT/AC).")
bam_to_elr_parser.add_argument("--sj_shift", dest='SJ_SHIFT', default=0, type=int, help="Shift up to this many bases to find a canonical splice junction")
bam_to_elr_parser.add_argument("--minlen_strict", dest='MINLEN_STRICT', default=18, type=int, help="Keep reads down to this length only if perfectly aligned.")
bam_to_elr_parser.add_argument("--minlen_loose", dest='MINLEN_LOOSE', default=25, type=int, help="Keep reads down to this length if they passed alignment parameters.")
bam_to_elr_parser.add_argument("--error_rate", dest='ERROR_RATE', default=0.10, type=float, help="Maximum allowed error rate (mismatches+indels) per exon.")
bam_to_elr_parser.add_argument("INPUT", type=str, default=None, help="Input BAM/SAM file")
bam_to_elr_parser.set_defaults(object='BAMtoELRconverter')

### bed_to_elr.py ###
bed_to_elr_parser = subparsers.add_parser('bed-to-elr',help="Converts a BED file to an End-Labeled Read (ELR) file.", description=ELRdesc, formatter_class=ArgumentDefaultsHelpFormatter)
bed_to_elr_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, required=True, help="Filepath to write ELR file.")
bed_to_elr_parser.add_argument("--chroms", dest='CHROMS', type=str, help="Filepath to a text file of chromosome names (1 per line).")
bed_to_elr_parser.add_argument("--header", dest='HEADER', type=str, default=None, help="Filepath to write ELR header.")
bed_to_elr_parser.add_argument("--source", dest='SOURCE', help="Source of BED lines.", default='BED', type=str)
bed_to_elr_parser.add_argument("-j", dest='JUNCTIONS', default=False, action='store_true', help="Gaps in the reads are all splice junctions.")
bed_to_elr_parser.add_argument("-s", dest='START', default=False, action='store_true', help="All read 5' ends are transcript start sites.")
bed_to_elr_parser.add_argument("-c", dest='CAPPED', default=False, action='store_true', help="All 5' ends are capped.")
bed_to_elr_parser.add_argument("-e", dest='END', default=False, action='store_true', help="All read 3' ends are transcript end sites.")
bed_to_elr_parser.add_argument("INPUT", type=str, help='Input BED file')
bed_to_elr_parser.set_defaults(object='BEDtoELRconverter')

### elr_combine.py ###
combine_parser = subparsers.add_parser('elr-combine',help="Makes one unified End-Labeled Read (ELR) file from multiple sorted files.", description=ELRdesc, formatter_class=ArgumentDefaultsHelpFormatter)
combine_parser.add_argument(dest='INPUT', help="Input sorted ELR files.", type=str, nargs='+')
combine_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, default=None, required=True, help="Filepath to write BED file.")
combine_parser.add_argument("--scale", dest='SCALE', type=float, nargs='+', default=None, help="Scaling factors (1 per input file)")
combine_parser.add_argument("--source", dest='SOURCE', default=None, type=str, help="Merge all input reads into a single source.")
combine_parser.add_argument("--temp", dest='TEMPDIR', type=str, default='_combinetmp', help="Prefix for temp files.")
combine_parser.set_defaults(object='ELRcombiner')


### elr_condense.py ###
condense_parser = subparsers.add_parser('condense',help="Partial assembly an end-labeled read (ELR) file. Outputs all loci (no filters) to a new sorted ELR.", formatter_class=ArgumentDefaultsHelpFormatter)
condense_parser.add_argument('-o','--output', dest='OUT', type=str, default=None, help="Destination file for assembly. File extension (bed, elr, gtf) determines output type.")
condense_parser.add_argument('--max_gap', dest='MAX_GAP', type=int, default=0, help="Largest gap size to tolerate (nucleotides).")
condense_parser.add_argument('--end_cluster', dest='END_CLUSTER', type=int, default=50, help="Largest distance between end-labeled reads to consider the same cluster (nucleotides).")
condense_parser.add_argument('--min_overhang', dest='MIN_OVERHANG', type=int, default=3, help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).")
condense_parser.add_argument('--min_cov', dest='MIN_COV', type=float, default=1, help="Minimum coverage filter to remove low-evidence transcript models.")
condense_parser.add_argument('--min_len', dest='MINLEN', type=int, default=20, help="Minimum output transcript length (nucleotides).")
condense_parser.add_argument('--min_intron_len', dest='MIN_INTRON_LEN', type=int, default=50, help="Minimum length for an intron (nucleotides)")
condense_parser.add_argument('--min_proportion', dest='MIN_PROPORTION', type=float, default=0.01, help="[float 0-1] Exclude ends, juctions, or transcripts that contribute < this proportion. (Used as a signal threshold)")
condense_parser.add_argument('--intron_filter', dest='INTRON_FILTER', type=float, default=0.15, help="[float 0-1] Retained introns must exceed this proportion the be considered.")
condense_parser.add_argument('--cap_bonus', dest='CAP_BONUS', type=float, default=5, help="[float] Signal multiplier for 5' reads with an inferred cap structure (uuG).")
condense_parser.add_argument('--cap_filter', dest='CAP_FILTER', type=float, default=.1, help="[float] Threshold percent uuG to count cluster as capped.")
condense_parser.add_argument("--starts", dest='STARTS', default=False, action='store_true', help="Sample is a Start Tag (5' end) file, e.g. CAGE")
condense_parser.add_argument("--ends", dest='ENDS', default=False, action='store_true', help="Sample is an End Tag (3' end) file, e.g. 3P-Seq")
condense_parser.add_argument("--sparse", dest='SPARSE', default=False, action='store_true', help="Sample is sparsely end-labeled, e.g. Smart-seq")
condense_parser.add_argument(dest='INPUT', type=str, help="Input single ELR filepath. MUST be coordinate-sorted.")
condense_parser.set_defaults(object='Condenser')


### elr_to_bed.py ###
elr_to_bed_parser = subparsers.add_parser('elr-to-bed',help="Converts an End-Labeled Read (ELR) file to BED12.", formatter_class=ArgumentDefaultsHelpFormatter)
elr_to_bed_parser.add_argument("-o", "--output", dest='OUTPUT', type=str, default=None, required=True, help="Filepath to write BED file.")
elr_to_bed_parser.add_argument("INPUT", help='Input ELR file')
elr_to_bed_parser.set_defaults(object='ELRtoBEDconverter')

### fasta.py ###
fasta_parser = subparsers.add_parser('fasta',help="Writes a transcript FASTA file for each input feature.", formatter_class=ArgumentDefaultsHelpFormatter)
fasta_parser.add_argument("-o", "--output", dest='OUT', type=str, default='bookend.fasta', help="Filepath to write feature FASTA file.")
fasta_parser.add_argument('--genome', dest='GENOME', required=True, help="(required) Path to genome FASTA file.", default=None)
fasta_parser.add_argument('--allow_unstranded', dest='UNSTRANDED', default=False, action='store_true', help="Allow unstranded transcripts to be written to output (forward strand).")
fasta_parser.add_argument('--orf', dest='ORF', default=False, action='store_true', help="Report longest open reading frame.")
fasta_parser.add_argument('--min_aa', dest='MIN_AA', type=int, default=10, help="Minimum ORF length (in amino acids) to report.")
fasta_parser.add_argument('--allow_partial_orf', dest='PARTIAL_ORF', default=False, action='store_true', help="Allow ORFs missing a start and/or stop codon.")
fasta_parser.add_argument('-w','--width', dest='WIDTH', type=int, default=80, help="Number of columns for the output FASTA.")
fasta_parser.add_argument('INPUT', type=str, help="[GFF3/GTF/ELR/BED] Path to feature file(s)", nargs='*')
fasta_parser.set_defaults(object='FastaWriter')


### fastq_end_label.py ###    
end_label_parser = subparsers.add_parser('label',help="Trims and labels RNA 5' and 3' ends in a FASTQ file", formatter_class=ArgumentDefaultsHelpFormatter)
end_label_parser.add_argument('-S', '--start', dest='START', type=str, default='AAGCAGTGGTATCAACGCAGAGTACATGGG', help="Template switching primer, marks the RNA 5' terminus.")
end_label_parser.add_argument( '-E', '--end', dest='END', type=str, default='ACGCAGAGTACTTTTTTTTTTTTTTTTTTTT+', help="Reverse transcription primer, marks (reverse complement of) the RNA 3' terminus.")
end_label_parser.add_argument('--umi', dest='UMI', type=str, default="", choices=['S', 'E'], help="One of the end labels (S, E) contains a Unique Molecular Identifier (UMI).")
end_label_parser.add_argument('--strand', dest='STRAND', type=str, default='unstranded', choices=['forward','reverse','unstranded'], help="Orientation of mate1 with respect to the RNA strand")
end_label_parser.add_argument('--out1', dest='OUT1', type=str, default='label.1.fq', help="Destination file for trimmed and labeled mate 1 FASTQ file.")
end_label_parser.add_argument('--out2', dest='OUT2', type=str, default='label.2.fq', help="Destination file for trimmed and labeled mate 2 FASTQ file.")
end_label_parser.add_argument('--single_out', dest='SINGLE_OUT', type=str, default='label.single.fq', help="Destination file for reads that are only informative single-end.")
end_label_parser.add_argument('--min_start', dest='MIN_START', type=int, default=7, help="Minimum number of nucleotides that must match to the start adapter sequence.")
end_label_parser.add_argument('--min_end', dest='MIN_END', type=int, default=9, help="Minimum number of nucleotides that must match to the end adapter sequence.")
end_label_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display each trimming result on stdout.")
end_label_parser.add_argument('--max_start', dest='MAX_START', type=int, default=120, help="Maximum number of nucleotides to trim for the start adapter sequence.")
end_label_parser.add_argument('--max_end', dest='MAX_END', type=int, default=120, help="Maximum number of nucleotides to trim for the end adapter sequence.")
end_label_parser.add_argument('--minlen', dest='MINLEN', type=int, default=18, help="Minimum sequence length to keep.")
end_label_parser.add_argument('--mismatch_rate', dest='MM_RATE', type=float, default=.06, help="Highest allow proportion of mismatches.")
end_label_parser.add_argument('--qualmask', dest='QUALMASK', type=float, default=16, help="Ignores any basecalls with phred score < this, treats base as N.")
end_label_parser.add_argument('--minqual', dest='MINQUAL', type=float, default=25, help="Suppresses any trimmed sequences with lower than this mean phred quality score.")
end_label_parser.add_argument('--discard_untrimmed', dest='DISCARD_UNTRIMMED', default=False, action='store_true', help="If no trimming occurred, do not write the read to output.")
end_label_parser.add_argument(dest='FASTQ', type=str, nargs='+', help="Input FASTQ file(s). 1 for single-end, 2 for paired-end")
end_label_parser.add_argument('--pseudomates', dest='PSEUDOMATES', default=False, action='store_true', help="Write single reads to --out1 with an artificial reverse complement mate pair in --out2 (overrides --single_out)")
end_label_parser.set_defaults(object='EndLabeler')

### gtf_merge.py ###
merge_parser = subparsers.add_parser('merge',help="Merges multiple assembly/annotation files into one consensus annotation.", formatter_class=ArgumentDefaultsHelpFormatter)
merge_parser.add_argument('-r', dest='REFERENCE', help="[GFF3/GTF] Path to reference annotation", type=str, default=None)
merge_parser.add_argument("-i", "--input", dest="INPUT", type=str, nargs='+', help="Input assembly GTF/GFF3/BED12 file(s)")
merge_parser.add_argument("-o", "--output", dest='OUT', type=str, default='merge.gtf', help="Filepath to write output class file")
merge_parser.add_argument('--name', dest='NAME', help="Prefix for output transcripts.", type=str, default='BOOKEND')
merge_parser.add_argument("--genome", dest='GENOME', default=None, type=str, help="Genome FASTA file. Required for checking splice motifs.")
merge_parser.add_argument('--keep_refs', dest='KEEP_REFS', default=False, action='store_true', help="Output all reference transcripts, even if unsupported.")
merge_parser.add_argument("--table", dest='TABLE', type=str, default='merge_summary.tsv', help="Filepath to write summary TSV file")
merge_parser.add_argument('--end_buffer', dest='END_BUFFER', type=int, default=100, help="Largest distance between ends to be considered the same (nucleotides).")
merge_parser.add_argument('--ref_parent', dest='GFF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Parent object")
merge_parser.add_argument('--min_terminal', dest='MIN_TERMINAL', type=int, default=40, help="Smallest allowed terminal exon.")
merge_parser.add_argument('--ref_child', dest='GFF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Child object")
merge_parser.add_argument('--parent_attr_gene', dest='PARENT_ATTR_GENE', type=str, default=None, help="Gene attribute in Parent objects")
merge_parser.add_argument('--child_attr_gene', dest='CHILD_ATTR_GENE', type=str,  default=None, help="Gene attribute in Child objects")
merge_parser.add_argument('--parent_attr_transcript', dest='PARENT_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Parent objects")
merge_parser.add_argument('--child_attr_transcript', dest='CHILD_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Child objects")
merge_parser.add_argument('--bed_gene_delim', dest='GENE_DELIM', type=str, default='.', help="(for BED12 references) String that splits gene name from transcript isoform.")
merge_parser.add_argument('--allow_unstranded', dest='UNSTRANDED', default=False, action='store_true', help="Allow unstranded transcripts to match overlapping transcripts.")
merge_parser.add_argument('--fusion_delim', dest='FUSION_DELIM', type=str, default='|', help="Character to use in fused transcript names.")
merge_parser.add_argument('--refname', dest='REFNAME', help="Name of the reference source", type=str, default=None)
merge_parser.add_argument('--min_len', dest='MIN_LEN', type=int, default=100, help="Minimum transcript length of merged assemblies.")
merge_parser.add_argument('--min_start', dest='MIN_S', type=float, default=0, help="Minimum number of start reads to accept as a start site.")
merge_parser.add_argument('--min_end', dest='MIN_E', type=float, default=0, help="Minimum number of end reads to accept as an end site")
merge_parser.add_argument('--attr_merge', dest='ATTR_MERGE', type=str, default='sum', choices=['sum','max'], help="Method for combining numeric attributes.")
merge_parser.add_argument('--rep_filter', dest='REP_FILTER', type=int, default=1, help="Number of times a transcript must be independently assembled (for multiple input files).")
merge_parser.add_argument('--tpm_filter', dest='TPM_FILTER', type=float, default=1, help="Ignore all transcripts with abundance lower than this.")
merge_parser.add_argument('--high_conf', dest='CONFIDENCE', type=float, default=3, help="[float] Multiplier for filters (reps, TPM) .")
merge_parser.add_argument('--cap_percent', dest='CAP_PERCENT', type=float, default=0, help="[float 0-1] Discard 5' end features with < this percent of cap-labeled reads (reads that contain upstream untemplated G).")
merge_parser.add_argument('--discard', dest='DISCARD', type=str, default=['ambiguous'], nargs='+', choices=['full_match', 'exon_match', 'isoform', 'fusion', 'fragment', 'intronic', 'antisense', 'intergenic', 'ambiguous'], help="Transcript class(es) to discard from merged annotation.")
merge_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each locus in stdout.")
merge_parser.set_defaults(object='AnnotationMerger')

### gtf_classify.py ###
classify_parser = subparsers.add_parser('classify',help="Classifies each transcript in an assembly against those in a reference annotation.", formatter_class=ArgumentDefaultsHelpFormatter)
classify_parser.add_argument("-i", "--input", dest="INPUT", type=str, nargs='+', help="Input assembly GTF/GFF3/BED12 file(s)")
classify_parser.add_argument("-o", "--output", dest='OUT', type=str, default='classify.tsv', help="Filepath to write output class file")
classify_parser.add_argument('--end_buffer', dest='END_BUFFER', type=int, default=100, help="Largest distance between ends to be considered the same (nucleotides).")
classify_parser.add_argument('-r', dest='REFERENCE', help="[GFF3/GTF] Path to reference annotation", type=str, default=None)
classify_parser.add_argument('--ref_parent', dest='GFF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Parent object")
classify_parser.add_argument('--ref_child', dest='GFF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Child object")
classify_parser.add_argument('--parent_attr_gene', dest='PARENT_ATTR_GENE', type=str, default=None, help="Gene attribute in Parent objects")
classify_parser.add_argument('--child_attr_gene', dest='CHILD_ATTR_GENE', type=str,  default=None, help="Gene attribute in Child objects")
classify_parser.add_argument('--parent_attr_transcript', dest='PARENT_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Parent objects")
classify_parser.add_argument('--child_attr_transcript', dest='CHILD_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Child objects")
classify_parser.add_argument('--pass_attr', dest='PASS_ATTR', default=False, action='store_true', help="Write attributes from reference transcripts instead of input transcripts.")
classify_parser.add_argument('--bed_gene_delim', dest='GENE_DELIM', type=str, default='.', help="(for BED12 references) String that splits gene name from transcript isoform.")
classify_parser.add_argument('--fusion_delim', dest='FUSION_DELIM', type=str, default='|', help="Character to use when joining gene names in a fused transcript.")
classify_parser.add_argument('--allow_unstranded', dest='UNSTRANDED', default=False, action='store_true', help="Allow unstranded transcripts to match overlapping transcripts.")
classify_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each locus in stdout.")
classify_parser.set_defaults(object='AssemblyClassifier')

### elr_quantify.py ###
quantify_parser = subparsers.add_parser('quantify',help="Quantifies the abundance of all annotated transcripts in each sample of the input ELR file(s).", formatter_class=ArgumentDefaultsHelpFormatter)
quantify_parser.add_argument("-i", "--input", dest="INPUT", type=str, nargs='+', help="Input assembly ELR file(s)")
quantify_parser.add_argument("-o", "--output", dest='OUT', type=str, default='quantify.tsv', help="Filepath to write output quantification file.")
quantify_parser.add_argument('-r', dest='REFERENCE', help="[GFF3/GTF/BED12] Path to reference annotation.", type=str, default=None, required=True)
quantify_parser.add_argument('--end_cluster', dest='END_CLUSTER', type=int, default=100, help="Largest distance between ends to be considered the same (nucleotides).")
quantify_parser.add_argument('--min_overhang', dest='MIN_OVERHANG', type=int, default=3, help="Smallest overhang to count for a read overlapping two exon fragments (number of nucleotides).")
quantify_parser.add_argument('--max_gap', dest='MAX_GAP', type=int, default=100, help="Largest gap allowed to process as a single chunk (nucleotides).")
quantify_parser.add_argument('--max_intron', dest='MAX_INTRON', type=int, default=100000, help="Ignore reads with introns longer than this.")
quantify_parser.add_argument('--ref_parent', dest='GFF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Parent object")
quantify_parser.add_argument('--ref_child', dest='GFF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in reference that define the Child object")
quantify_parser.add_argument('--parent_attr_gene', dest='PARENT_ATTR_GENE', type=str, default=None, help="Gene attribute in Parent objects")
quantify_parser.add_argument('--child_attr_gene', dest='CHILD_ATTR_GENE', type=str,  default=None, help="Gene attribute in Child objects")
quantify_parser.add_argument('--parent_attr_transcript', dest='PARENT_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Parent objects")
quantify_parser.add_argument('--child_attr_transcript', dest='CHILD_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Child objects")
quantify_parser.add_argument('--bed_gene_delim', dest='GENE_DELIM', type=str, default='.', help="(for BED12 references) String that splits gene name from transcript isoform.")
quantify_parser.add_argument('--verbose', dest='VERBOSE', default=False, action='store_true', help="Display a verbose summary of each locus in stdout.")
quantify_parser.set_defaults(object='Quantifier')


### sam_sj_out.py ###
sam_sj_parser = subparsers.add_parser('sam-to-sj',help="Generates a splice junction file (SJ.out.tab) from SAM.", formatter_class=ArgumentDefaultsHelpFormatter)
sam_sj_parser.add_argument("-F", "--fasta", dest='FASTA', help="Genome FASTA file", default=None, type=str, required=True)
sam_sj_parser.add_argument("--format", dest='FORMAT', help="Output file format", default='star', type=str, choices=['bed','star'])
sam_sj_parser.add_argument("--filter", dest='FILTER', help="Remove noncanonical splice junctions from the output", default=False, action='store_true')
sam_sj_parser.add_argument("INPUT", type=str, help="Input SAM file")
sam_sj_parser.set_defaults(object='SAMtoSJconverter')

### sj_merge.py ###
sj_merge_parser = subparsers.add_parser('sj-merge', help="Combines multiple SJ.out.tab or SJ.bed files.", formatter_class=ArgumentDefaultsHelpFormatter)
sj_merge_parser.add_argument("-o", "--output", dest='OUT', type=str, default='sj_merge.out.tab', help="Filepath to write merged file.")
sj_merge_parser.add_argument("--format", dest='FORMAT', help="Output file format", default='star', type=str, choices=['bed','star'])
sj_merge_parser.add_argument("--min_unique", dest='MIN_UNIQUE', help="Filter SJs with fewer unique reads.", default=0, type=int)
sj_merge_parser.add_argument("--min_reps", dest='MIN_REPS', help="Filter SJs detected in fewer than this many files.", default=1, type=int)
sj_merge_parser.add_argument("--new", dest='NEW', help="Keep only SJs not present in the reference SJDB", default=False, action='store_true')
sj_merge_parser.add_argument("INPUT", nargs='+', default=[])
sj_merge_parser.set_defaults(object='SJmerger')

### sj_to_bed.py ###
sj_to_bed_parser = subparsers.add_parser('sj-to-bed', help="Converts SJ.out.tab file to an SJ.bed file.", formatter_class=ArgumentDefaultsHelpFormatter)
sj_to_bed_parser.add_argument("INPUT", type=str, help="Input SJ.out.tab file")
sj_to_bed_parser.add_argument("-o", "--output", dest='OUT', type=str, default='SJ.bed', help="Filepath to write SJ.bed file.")
sj_to_bed_parser.set_defaults(object='SJtoBEDconverter')

### elr_sort.py ###
elr_sort_parser = subparsers.add_parser('elr-sort',help="Sorts an End-Labeled Read (ELR) file.", formatter_class=ArgumentDefaultsHelpFormatter)
elr_sort_parser.add_argument("-o", "--output", dest='OUT', help="Output file path (default: stdout)", default='stdout')
elr_sort_parser.add_argument("-f" ,"--force", dest='FORCE', help="Force overwrite of --output file if it exists.", default=False, action='store_true')
elr_sort_parser.add_argument("INPUT", type=str, help="Input ELR file")
elr_sort_parser.set_defaults(object='ELRsorter')

### elr_subset.py ###
elr_subset_parser = subparsers.add_parser('elr-subset',help="Writes a subsetted region of an ELR file.", formatter_class=ArgumentDefaultsHelpFormatter)
elr_subset_parser.add_argument("-o", "--output", dest='OUT', help="Output file path (default: stdout)", default='stdout')
elr_subset_parser.add_argument("-f" ,"--force", dest='FORCE', help="Force overwrite of --output file if it exists.", default=False, action='store_true')
elr_subset_parser.add_argument("-r" ,"--region", dest='REGION', help="[chrom:start-end] Region to write to output", type=str, required=True)
elr_subset_parser.add_argument("INPUT", type=str, help="Input ELR file")
elr_subset_parser.set_defaults(object='ELRsubsetter')


### gtf_to_bed.py ###
gtf_to_bed_parser = subparsers.add_parser('gtf-to-bed',help="Converts a GTF/GFF3 annotation file to BED12.", formatter_class=ArgumentDefaultsHelpFormatter)
gtf_to_bed_parser.add_argument("INPUT", type=str, help="Input GTF/GFF3 file")
gtf_to_bed_parser.add_argument("-o", "--output", dest='OUT', help="Output file path (default: stdout)", default='stdout')
gtf_to_bed_parser.add_argument("-f" ,"--force", dest='FORCE', help="Force overwrite of --output file if it exists.", default=False, action='store_true')
gtf_to_bed_parser.add_argument("--source", dest='SOURCE', help="Source of GTF/GFF3 file (default: file name)", default=None, type=str)
gtf_to_bed_parser.add_argument("--name", dest='NAME_ATTR', help="Attribute to pass to the name column (default: transcript_id)", default='transcript_id', type=str)
gtf_to_bed_parser.add_argument("--score", dest='SCORE', help="Attribute to pass to score column (default: None)", default=None, type=str)
gtf_to_bed_parser.add_argument('--gtf_parent', dest='GTF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in GTF files for Parent object (default: transcript)")
gtf_to_bed_parser.add_argument('--gtf_child', dest='GTF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in GTF files for Child object (default: exon)")
gtf_to_bed_parser.add_argument('--gff_parent', dest='GFF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in GFF3 files for Parent object (default: mRNA, transcript)")
gtf_to_bed_parser.add_argument('--gff_child', dest='GFF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in GFF3 files for Child object (default: exon)")
gtf_to_bed_parser.add_argument('--parent_attr_gene', dest='PARENT_ATTR_GENE', type=str, default=None, help="Gene attribute in Parent objects (default: gene)")
gtf_to_bed_parser.add_argument('--child_attr_gene', dest='CHILD_ATTR_GENE', type=str, default=None, help="Gene attribute in Child objects (default: gene)")
gtf_to_bed_parser.add_argument('--parent_attr_transcript', dest='PARENT_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Parent objects (default: transcript_id)")
gtf_to_bed_parser.add_argument('--child_attr_transcript', dest='CHILD_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Child objects (default: transcript_id)")
gtf_to_bed_parser.add_argument('--color_key', dest='COLOR_KEY', type=str, default=None, help="Attribute name to lookup transcript type")
gtf_to_bed_parser.add_argument('--color_code', dest='COLOR_CODE', type=str, default=None, help="Tab-separated file of transcript type R,G,B colors (e.g. 'type\t255,255,255')")
gtf_to_bed_parser.add_argument("--extend", dest='EXTEND', help="Extend annotation to ends of Start and End clusters, if they exist.", default=False, action='store_true')
gtf_to_bed_parser.set_defaults(object='GTFconverter')

### gtf_ends.py ###
gtf_ends_parser = subparsers.add_parser('gtf-ends',help="Writes a BED file containing the set of end features in an annotation (GTF/GFF3/BED12).", formatter_class=ArgumentDefaultsHelpFormatter)
gtf_ends_parser.add_argument("INPUT", type=str, help="Input GTF/GFF3 file")
gtf_ends_parser.add_argument("-o", "--output", dest='OUT', help="Output file path (default: stdout)", default='stdout')
gtf_ends_parser.add_argument("-f" ,"--force", dest='FORCE', help="Force overwrite of --output file if it exists.", default=False, action='store_true')
gtf_ends_parser.add_argument("--no_merge", dest='NO_MERGE', help="Do not collapse overlapping ends into a single BED feature.", default=False, action='store_true')
gtf_ends_parser.add_argument('-t','--type', dest='TYPE', choices=['all','S', 'E'], type=str, default='all', help="End type: Starts (S), Ends (E), or all")
gtf_ends_parser.add_argument("--extend", dest='EXTEND', type=int, help="Extend a fixed distance up and down from the peak position.", default=None)
gtf_ends_parser.add_argument('--min_terminal', dest='MIN_TERMINAL', type=int, default=40, help="Smallest allowed terminal exon.")
gtf_ends_parser.add_argument("--name", dest='NAME_ATTR', help="Attribute to pass to the name column (default: transcript_id)", default='transcript_id', type=str)
gtf_ends_parser.add_argument("--score", dest='SCORE', help="Attribute to pass to score column (default: None)", default=None, type=str)
gtf_ends_parser.add_argument('--gtf_parent', dest='GTF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in GTF files for Parent object (default: transcript)")
gtf_ends_parser.add_argument('--gtf_child', dest='GTF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in GTF files for Child object (default: exon)")
gtf_ends_parser.add_argument('--gff_parent', dest='GFF_PARENT', type=str, nargs='+', default=None, help="Line type(s) in GFF3 files for Parent object (default: mRNA, transcript)")
gtf_ends_parser.add_argument('--gff_child', dest='GFF_CHILD', type=str, nargs='+', default=None, help="Line type(s) in GFF3 files for Child object (default: exon)")
gtf_ends_parser.add_argument('--parent_attr_gene', dest='PARENT_ATTR_GENE', type=str, default=None, help="Gene attribute in Parent objects (default: gene)")
gtf_ends_parser.add_argument('--child_attr_gene', dest='CHILD_ATTR_GENE', type=str,  default=None, help="Gene attribute in Child objects (default: gene)")
gtf_ends_parser.add_argument('--parent_attr_transcript', dest='PARENT_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Parent objects (default: transcript_id)")
gtf_ends_parser.add_argument('--child_attr_transcript', dest='CHILD_ATTR_TRANSCRIPT', nargs='+', type=str, default=None, help="Transcript attribute(s) in Child objects (default: transcript_id)")
gtf_ends_parser.set_defaults(object='GTFendwriter')


### elr_simulate.py ###
# simulate_parser = subparsers.add_parser('simulate',help="Simulates an ELR file from a ground truth of transcripts and abundances.")
# simulate_parser.add_argument("-o", "--output", dest='OUT', help="Output file path (default: stdout)", default='stdout')
# simulate_parser.add_argument('-r', '--reference', dest='REFERENCE', help="Reference annotation file (ELR/BED12/GTF/GFF3).", default=None, type=str)
# simulate_parser.add_argument('-a', '--abundance', dest='ABUNDANCE', help="Two-column table of [transcript ID]\t[relative abundance].", default=None, type=str)
# simulate_parser.add_argument("--count", dest='READ_COUNT', type=int, default=1000000, help="Number of reads to simulate.")
# simulate_parser.add_argument("--read_length", dest='READ_LENGTH', type=int, default=50, help="Sequencing read length for simulated reads.")
# simulate_parser.add_argument("--min_length", dest='MIN_LENGTH', type=int, default=20, help="Smallest length of genome-matching read to keep.")
# simulate_parser.add_argument("--adapter_5p", dest='ADAPTER_5P', type=int, default=30, help="Length of adapter on the 5' end of cDNA.")
# simulate_parser.add_argument("--adapter_3p", dest='ADAPTER_3P', type=int, default=55, help="Length of adapters on the 3' end of cDNA.")
# simulate_parser.add_argument("--paired", dest='PAIRED', default=False, action="store_true", help="If True, sequencing is paired-end.")
# simulate_parser.add_argument("--fragment_mean", dest='FRAGMENT_MEAN', type=float, default=200, help="Mean length of RNA fragments.")
# simulate_parser.add_argument("--fragment_sd", dest='FRAGMENT_SD', type=float, default=50, help="Standard deviation of RNA fragment length (gaussian).")
# simulate_parser.add_argument("--var_5p", dest='VAR_5P', type=float, default=0, help="Variation of 5P end positions (laplace scale parameter).")
# simulate_parser.add_argument("--var_3p", dest='VAR_3P', type=float, default=0, help="Variation of 3P end positions (laplace scale parameter).")
# simulate_parser.add_argument("--percent_intact", dest='PERCENT_INTACT', type=float, default=100, help="(0-100] Percent of transcripts that are full-length.")
# simulate_parser.add_argument("--percent_sense", dest='PERCENT_SENSE', type=float, default=50, help="[0-100] Percent of reads that are in sense to the RNA strand.")
# simulate_parser.add_argument("--label_noise", dest='LABEL_NOISE', type=float, default=0, help="[0-100] Percent of reads falsely labeled ends.")
# simulate_parser.add_argument("--seed", dest='SEED', type=int, default=0, help="Seed for random number generator.")
# simulate_parser.add_argument("--bed", dest='BED', default=False, action="store_true", help="If True, output in BED12 format.")
# simulate_parser.set_defaults(object='ELRsimulator')
