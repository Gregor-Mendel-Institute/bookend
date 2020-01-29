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
        "Converts 'end labeled read' (ELR) files to BED.\n"
        "Three additional columns are added to the traditional 12-column format:\n"
        "\n"
        "  weight: Read count (allows partial counts for multimappers)\n"
        "  sample: Source of the ELR entry\n"
        "  end label: String describing the ends of each block as a start, end, junction, or none\n"
    )
    # initilize argumentparser to read in commands
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-O", "--output", dest='OUTPUT', type=str, default='stdout',
        help="Filepath to write ELR file."
    )
    parser.add_argument(
        "FILENAME", nargs='?'
    )
    args = parser.parse_args()
    
    s={'+':'plus','-':'minus','.':'ns'}
    if args.FILENAME:
        if args.FILENAME.split('.')[-1].lower() not in ['elr']:
            print("\nERROR: input file must be ELR format.")
            parser.print_help()
            sys.exit(1)
        input_file = open(args.FILENAME)
    elif not sys.stdin.isatty():
        input_file = sys.stdin
    else:
        print("\nERROR: requires ELR file as input.")
        parser.print_help()
        sys.exit(1)

    output_file = 'stdout'
    if args.OUTPUT != 'stdout':
        output_file = open(args.OUTPUT, 'w')
    
    #TODO: Allow dataset to be initialized with a separate header file
    dataset = ru.RNAseqDataset()
    
    for elr_line in input_file:
        if elr_line[0] == '#':
            header_line = elr_line.rstrip().split(' ')
            if header_line[0] == '#S':
                dataset.add_source(header_line[-1])
            if header_line[0] == '#C':
                dataset.add_chrom(header_line[-1])
            
            output_lines([elr_line.rstrip()], output_file)
            continue
        
        dataset.add_read_from_ELR(elr_line)
        while len(dataset.read_list) > 0:
            bed_line = dataset.pop_read('bed')
            output_lines([bed_line], output_file)
