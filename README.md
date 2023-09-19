# /| bookend |\\  
### End-guided transcriptome assembly.  
Bookend is a comprehensive framework for end-guided assembly of short-read, long-read, and end-capture RNA-seq data.
Please see the [User Guide](Bookend_User_Guide.pdf) for a full description of the subcommands and arguments.
The lastest developments can be found in the [Bookend GitHub repository](https://github.com/Gregor-Mendel-Institute/bookend).

Take a look at the [Bookend publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02700-3) for more details about its usage and applications.

## Installation  
Bookend can be installed through the Python Package Index (PyPI) on UNIX systems with Python 3.6+ using the command:
```
pip install bookend-rna --upgrade
```

If installing from the GitHub source code, perform the following steps:
```
git clone https://github.com/Gregor-Mendel-Institute/bookend
cd bookend
python3 setup.py install
```

Once installed, all utilities can be accessed on the command as bookend subcommands:  
  
    usage: bookend [subcommand] [options] [input file(s)]
    Subcommands (use -h/--help for more info):
    
        label    (Label 5' and 3' ends in a FASTQ file)
        elr      (Convert a BAM/SAM file to the end-labeled read format)
        assemble (Assemble transcripts from aligned end-labeled reads)
        condense (Partial assembly that leaves keeps all fragments; use for meta-assembly)
        classify (Compare an assembly to the transcripts of a reference annotation)
        merge    (Combine GTF files into a unified annotation)
        bedgraph (Write a coverage Bedgraph file of end-labeled reads)
        fasta    (Write a transcript FASTA file from an annotation and genome)
    
        --end-labeled read (ELR) operations--
        elr-sort
        elr-subset
        elr-combine
        
        --file conversion--
        gtf-to-bed
        gtf-ends
        bed-to-elr
        elr-to-bed
        sam-to-sj
        sj-to-bed
        sj-merge  
