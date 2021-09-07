# /| bookend |\\  
### End-guided transcriptome assembly.  
Bookend is a comprehensive framework for end-guided assembly of short-read, long-read, and end-capture RNA-seq data.
Please see the [User Guide](Bookend_User_Guide.pdf) for a full description of the subcommands and arguments.

## Installation  
Bookend can be installed through the Python Package Index (PyPI) on UNIX systems with Python 3.6+ using the command  
    pip install bookend-rna

Once installed, all utilities can be accessed on the command as bookend subcommands:  

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
  
