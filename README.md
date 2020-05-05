# /| bookend |\\  
### End-guided transcriptome assembly.

usage: bookend [subcommand] [options] [input file(s)]
Subcommands (use -h/--help for more info):

    label    (Label 5' and 3' ends in a FASTQ file)
    assemble (Assemble transcripts from aligned end-labeled reads)
    merge    (Merge assembled GTFs with a reference annotation)

    --end-labeled read (ELR) operations--
    make-elr
    sort-elr
    combine-elr

    --file conversion--
    bed-to-elr
    elr-to-bed
    gtf-to-bed
    sam-to-sj
    sj-to-bed
    sj-merge

    --setup/indexing--
    index-fasta
    softbridge-fasta
  
    
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