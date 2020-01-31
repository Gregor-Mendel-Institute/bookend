import os
import re
import sys
import argparse
import time
from cython_utils import _fasta_utils as fu

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    "-F", "--fasta", dest='FASTA',
    help="Genome FASTA file",
    default=None, type=str, required=True
)
parser.add_argument(
    "-O", "--output", dest='OUTPUT', type=str, default='stdout',
    help="Filepath to write BED12 file."
)
parser.add_argument(
    "--maxlen", dest='MAXLEN',
    help="Maximum tolerated length for a bridge over a softmasked region.",
    default=1000, type=int
)
args = parser.parse_args()

################################
# ENVIRONMENT SETUP: FUNCTIONS #
################################
    
def output_lines(lines, output):
    """Takes a list of bed lines and writes
    them to the output stream.
    """
    if output == 'stdout':
        for output_string in lines:
            print(output_string)
    else:
        for output_string in lines:
            output.write('{}\n'.format(output_string))

def generate_softbridges(genome_dict, maxlen):
    """From a genome dict, writes a generator
    that yields one BED12 line for each start/stop
    of a softmasked region of the FASTA file,
    demarcated by lowercase letters.
    """
    for chrom in sorted(list(genome_dict.keys())):
        soft_toggle = False
        current_pos = 0
        start_pos = 0
        end_pos = 0
        for lowercase in [i.islower() for i in genome_dict[chrom]]:
            if lowercase:
                if not soft_toggle:
                    # Start a softbridge
                    soft_toggle = True
                    start_pos = current_pos
            else:
                if soft_toggle:
                    # End a softbridge
                    soft_toggle = False
                    end_pos = current_pos
                    out_string = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        chrom,start_pos,end_pos,
                        'UU',0,'.',0,0,'204,204,180',1,end_pos-start_pos,0,
                        0.01,'softbridge','.'
                    )
                    yield out_string
            
            current_pos += 1

###########
# EXECUTE #
###########

if __name__ == '__main__':
    if args.FASTA:
        genome = fu.import_genome(args.FASTA)
    else:
        print('ERROR: Requires a reference genome. Use the -F/--fasta argument.')
        sys.exit(1)

    if args.OUTPUT != 'stdout':
        output_file = open(args.OUTPUT, 'w')
    else:
        output_file = 'stdout'
    
    output_lines(generate_softbridges(genome, args.MAXLEN), output_file)
    
    if output_file != 'stdout':
        output_file.close()
