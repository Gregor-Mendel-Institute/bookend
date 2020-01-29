import sys
import argparse
import _rnaseq_utils as ru
from argparse import RawTextHelpFormatter

def output_lines(lines, output):
    """Takes a list of bed lines and writes
    them to the output stream.
    """
    if output == 'stdout':
        for output_string in lines:
            print(output_string)
    else:
        for output_string in lines:
            o = output_string.rstrip()
            output.write('{}\n'.format(o))

if __name__ == "__main__":
    ###################
    # INPUT ARGUMENTS #
    ###################
    desc = (
        "Converts BED files to a 'end labeled read' (ELR) format.\n"
        "ELR contains a two-component header: #G (genome) and #N (names of samples).\n"
        "Each line is a read or read stack with six columns:\n"
        "  chromosome  position  strand  ELCIGAR  sample  weight\n"
        "\n"
        "  chromosome: Chromosome number, as indexed by the #G header\n"
        "  position: 0-indexed position in chromosome\n"
        "  strand: Inferred RNA strand; +, -, or . (unknown)\n"
        "  ELCIGAR: String describing the mapped read (full description below)\n"
        "  sample: Sample number, as indexed by the #N header\n"
        "  weight: Read count (allows partial counts for multimappers)\n"
    )

    epilog = (
        "ELCIGAR strings are Character|Number strings and one trailing character\n"
        "([CN]xC), where C is a label and N is a numeric length on the genome.\n"
        "Each C labels an end or gap in the alignment as one of the following:\n"
        "  S: start\n"
        "  s: start (low confidence)\n"
        "  C: start (capped)\n"
        "  E: end (polyA tail)\n"
        "  e: end (low confidence)\n"
        "  D: splice junction donor\n"
        "  A: splice junction acceptor\n"
        "  .: unspecified gap or end\n"
        "\n"
        "For example, a paired-end 50 read of a 185-nt fragment would be\n"
        "  x50x85x50x\n"
        "A full 3-exon transcript could be described as:\n"
        "  S256D800A128D800A512E\n"
        "where the 3 exons are 256, 128, and 512 nucleotides,\n"
        "and the 2 introns are both 800 nucleotides.\n"
        "\n"
    )

    # initilize argumentparser to read in commands
    parser = argparse.ArgumentParser(
        description=desc,
        epilog=epilog,
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "--source", dest='SOURCE',
        help="Source of BED lines.",
        default=None, type=str
    )
    parser.add_argument(
        "-j", dest='JUNCTIONS', default=False, action='store_true',
        help="Gaps in the reads are splice junctions."
    )
    parser.add_argument(
        "-s", dest='START', default=False, action='store_true',
        help="Read 5' ends are transcript start sites."
    )
    parser.add_argument(
        "-c", dest='CAPPED', default=False, action='store_true',
        help="5' end data is capped."
    )
    parser.add_argument(
        "-e", dest='END', default=False, action='store_true',
        help="Read 3' ends are transcript end sites."
    )
    parser.add_argument(
        "--bed", dest='BED_OUT', default=False, action='store_true',
        help="Output a 15-column end labeled BED file."
    )
    parser.add_argument(
        "-O", "--output", dest='OUTPUT', type=str, default='stdout',
        help="Filepath to write ELR file."
    )
    parser.add_argument(
        "--header", dest='HEADER', type=str, default='stdout',
        help="Filepath to write ELR header."
    )
    parser.add_argument(
        "FILENAME", nargs='?'
    )
    args = parser.parse_args()
    
    s={'+':'plus','-':'minus','.':'ns'}
    if args.FILENAME:
        if args.FILENAME.split('.')[-1].lower() not in ['bed','bed']:
            print("\nERROR: input file must be BED format.")
            parser.print_help()
            sys.exit(1)
        bed_in = open(args.FILENAME)
    elif not sys.stdin.isatty():
        bed_in = sys.stdin
    else:
        print("\nERROR: requires BED file as input.")
        parser.print_help()
        sys.exit(1)

    output_file = 'stdout'
    if args.OUTPUT != 'stdout':
        output_file = open(args.OUTPUT, 'w')
    
    header_file = 'stdout'
    if args.HEADER != 'stdout':
        header_file = open(args.HEADER, 'w')
    
    if args.SOURCE:
        source_array = [args.SOURCE]
    else:
        source_array = None
    
    dataset = ru.RNAseqDataset(source_array=source_array)
    current_source_index = dataset.source_index
    current_chrom_index = dataset.chrom_index
    
    if current_chrom_index > 0:
        output_lines(['#C {} {}'.format(i,c) for i,c in enumerate(dataset.chrom_array)], header_file)                      
                      
    if current_source_index > 0:
        output_lines(['#S {} {}'.format(i,c) for i,c in enumerate(dataset.source_array)], header_file)
    
    for bed_line in bed_in:
        if bed_line[0] == '#': # A header is being passed from the ELR file
            header_line = bed_line.rstrip().split(' ')
            if header_line[0] == '#S':
                dataset.add_source(header_line[-1])
            if header_line[0] == '#C':
                dataset.add_chrom(header_line[-1])
            
            output_lines([bed_line.rstrip()], output_file)
            continue
        
        current_chrom_index = dataset.chrom_index
        current_source_index = dataset.source_index
        dataset.add_read_from_BED(bed_line, source_string=args.SOURCE, s_tag=args.START, e_tag=args.END, capped=args.CAPPED, gaps_are_junctions=args.JUNCTIONS)
        if dataset.chrom_index > current_chrom_index: # A new chromosome was encountered
            output_lines(['#C {} {}'.format(len(dataset.chrom_array)-1, dataset.chrom_array[-1])], header_file)
            current_chrom_index = dataset.chrom_index
        
        if dataset.source_index > current_source_index: # A new source was encountered
            output_lines(['#S {} {}'.format(len(dataset.source_array)-1, dataset.source_array[-1])], header_file)
            current_source_index = dataset.source_index
        
        while len(dataset.read_list) > 0:
            if args.BED_OUT:
                out_line = dataset.pop_read('bed')
            else:
                out_line = dataset.pop_read()
            
            output_lines([out_line], output_file)



