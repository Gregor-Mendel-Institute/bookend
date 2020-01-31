import os
import re
import sys
import argparse
import time

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    'FASTA', type=str,
    help="Path to FASTA file."
)
parser.add_argument(
    '-S','--split', dest='SPLIT_ON', type=str,
    help="Character string to split header name from description.",
    default=' '
)
args = parser.parse_args()

################################
# ENVIRONMENT SETUP: FUNCTIONS #
################################
    
chrom='none'
FASTA_sequence = ''
chromosomes={}
trimmed_length = 'none'
genome_file=open(args.FASTA)
while genome_file:
    line = genome_file.readline()
    if line:
        l=line.rstrip()
        if l[0]=='>':
            if chrom!='none':
                FASTA_sequence = ''.join(x)
                print('{}\t{}\t{}\t{}\t{}'.format(
                    chrom,
                    len(FASTA_sequence),
                    start_position,
                    trimmed_length,
                    untrimmed_length
                ))
            chrom=l[1:len(l)].split(args.SPLIT_ON)[0]
            start_position = genome_file.tell()
            trimmed_length = 'none'
            length_mismatch = False
            x=[]
            continue
        else:
            if length_mismatch:
                stop("ERROR: [{}] line length mismatch".format(genome_file.tell()))
                
            if trimmed_length == 'none':
                trimmed_length = len(l)
                untrimmed_length = len(line)
            elif trimmed_length != len(l):
                length_mismatch = True
            
        x.append(l)
    else:
        FASTA_sequence = ''.join(x)
        print('{}\t{}\t{}\t{}\t{}'.format(
            chrom,
            len(FASTA_sequence),
            start_position,
            trimmed_length,
            untrimmed_length
        ))
        genome_file.close()
        genome_file = False
        
