import os
import sys
import argparse
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
    "--format", dest='FORMAT',
    help="Output file format",
    default='star', type=str, choices=['bed','star']
)
parser.add_argument(
    "--filter", dest='FILTER',
    help="Remove noncanonical splice junctions from the output",
    default=False, action='store_true'
)
parser.add_argument(
    "FILENAME", nargs='?'
)
args = parser.parse_args()

# A set of functions for updating an RNASeqRead object based on the genome FASTA file
def get_flank(chrom, pos, strand, label_type, label_len):
    """Gets flanking region to a genomic position to compare to an end label
    """
    flank = 'X'*label_len
    if label_type == 'S':
        if strand == '+':
            flank = genome[chrom][pos-label_len:pos]
        else:
            flank = fu.rc(genome[chrom][pos+1:pos+1+label_len])
    if label_type == 'E':
        if strand == '+':
            flank = genome[chrom][pos+1:pos+1+label_len]
        else:
            flank = fu.rc(genome[chrom][pos-label_len:pos])
    
    return flank

junction_types = {'GTAG':1, 'CTAC':2, 'GCAG':3, 'CTGC':4, 'ATAC':5, 'GTAT':6}    
def get_junction_strand(chrom,left,right):
    """ Returns '+', '-', or '.' for a left/right
    pair of splice junction positions based
    on the flanking genomic sequence """
    flanking_sequence = get_flank(chrom,left-1,'+','E',2) + get_flank(chrom,right,'+','S',2)
    flanking_sequence = flanking_sequence.upper()
    jtype = junction_types.get(flanking_sequence,0)
    return jtype

def parse_SAM_CIGAR(pos,cigar):
    """Converts the pos+CIGAR string of a SAM file to 
    an array of 0-indexed open (left,right) blocks as exons and introns
    (description from http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/)
    CIGAR operators:
    D    Deletion; present in the reference but not in the read
    H    Hard Clipping; not present in the read
    I    Insertion; present in the read  but not in the reference
    M    Match; can be either an alignment match or mismatch
    N    Skipped region; a region is not present in the read
    P    Padding; padded area in the read and not in the reference
    S    Soft Clipping;  the clipped nucleotides are present in the read
    """
    exons = []
    introns = []
    head = 0
    tail = 0
    current_pos = pos - 1
    cigar_breaks = [i for i,char in enumerate(cigar) if not char.isdigit()]
    l = 0
    first_element = True
    extend = False
    for b in cigar_breaks:
        operator, number, l = cigar[b], int(cigar[l:b]), b+1 # Get next CIGAR element
        if operator == 'S': # Softclipped ranges are stored in either head or tail
            if first_element:
                head = number
            else:
                tail = number
        
        first_element = False
        if operator == 'M':
            if extend: # Extends the previous blocks right side
                exons[-1] = (exons[-1][0], current_pos+number)
            else:
                exons += [(current_pos,current_pos+number)]
            
            current_pos += number
            extend = False
        elif operator == 'N':
            leftside = current_pos
            current_pos += number
            rightside = current_pos
            introns += [(leftside,rightside)]
        elif operator == 'D':
            current_pos += number
            extend = True
        elif operator == 'I':
            extend = True
    
    return exons, introns, head, tail

def junctions_from_sam(input_lines, sjdb):
    """Extract all splice junctions from line(s) of a SAM file"""
    ID = None
    if len(input_lines) > 1:
        multiline = True
    else:
        multiline = False
    
    for line in input_lines:
        l = line.strip().split('\t')
        line_ID = l[0]
        if not ID:
            ID = line_ID
        
        assert line_ID == ID, 'ERROR: nonmatching IDs in input_lines:\n{}\t{}'.format(ID,line_ID)
        
        SAMflags   = bin(int(l[1]))[2:]
        SAMflags   = '0'*(12-len(SAMflags))+SAMflags
        # Interpret the binary SAM flags
        is_paired      = bool(int(SAMflags[-1]))
        unmapped       = bool(int(SAMflags[-3]))
        read_reverse   = bool(int(SAMflags[-5]))
        first_in_pair  = bool(int(SAMflags[-7]))
        secondary      = bool(int(SAMflags[-9]))
        supplementary  = bool(int(SAMflags[-12]))
        
        chrom      = l[2]
        pos        = int(l[3])
        cigar      = l[5]
        seq        = l[9]
        
        if unmapped or supplementary:
            continue
        
        if multiline and not is_paired:
            multimapper = True
        else:
            if secondary:
                multimapper = True
            else:
                multimapper = False
        
        if read_reverse:
            strand = '-'
        else:
            strand = '+'
        
        if first_in_pair or not is_paired:
            mate = 1
        else: # Read is mate2 of a pair
            mate = 2
            if strand == '+': # Flip the strand to reflect the pair's strand, not the individual read
                strand = '-'
            else:
                strand = '+'
        
        # Populate sjdb with the current mapping
        if chrom not in sjdb:
            sjdb[chrom] = {}
        
        if strand == '+':
            s = 1
        elif strand == '-':
            s = 2
        else:
            s = 0
        
        if multimapper:
            unique, multi = 0, 1
        else:
            unique, multi = 1, 0
        
        exons, introns, head, tail = parse_SAM_CIGAR(pos, cigar)
        left_overhangs = [b-a for a,b in exons[:-1]]
        right_overhangs = [b-a for a,b in exons[1:]]
        for intron, lo, ro in zip(introns,left_overhangs,right_overhangs):
            motif = get_junction_strand(chrom,intron[0],intron[1])
            left, right = intron[0]+1, intron[1] # Convert to 1-indexed closed coords
            anno = 0
            overhang = min(lo,ro)
            
            feature_hash = '{}:{}:{}:{}:{}'.format(left,right,s,motif,anno)
            x = sjdb[chrom].get(feature_hash,(0,0,0))
            y = (int(unique),int(multi),int(overhang))
            sjdb[chrom][feature_hash] = (x[0]+y[0],x[1]+y[1],max([x[2],y[2]]))

######################
# PROCESS INPUT FILE # 
######################
if __name__ == '__main__':
    genome = fu.import_genome(args.FASTA)
    sj_dict = {}
    
    if args.FILENAME:
        if args.FILENAME.split('.')[-1].lower() != 'sam':
            print("\nERROR: input file must be SAM format.")
            parser.print_help()
            sys.exit(1)
        sam_in = open(args.FILENAME)
    elif not sys.stdin.isatty():
        sam_in = sys.stdin
    else:
        print("\nERROR: requires SAM file as input.")
        parser.print_help()
        sys.exit(1)

    
    chromosomes = {}
    header_written = False
    current_ID    = None
    current_read  = None
    current_lines = []
    for line in sam_in:
        if line[0] == '@':
            continue
        
        l = line.rstrip().split('\t')
        ID = l[0]
        if ID == current_ID or current_ID is None:
            current_lines.append(line)
            current_ID = ID
        else:
            # A new read ID was encountered. Perform actions on current_read
            # and begin a new current_read with the new ID.
            junctions_from_sam(current_lines,sj_dict)
            current_lines = [line]
            current_ID = ID

    junctions_from_sam(current_lines,sj_dict)
    
    # After the input is read all the way through, print sj_dict
    for chrom in sorted(list(sj_dict.keys())):
        if args.FORMAT == 'star':
            for item in sorted([tuple([int(i) for i in coords.split(':')])+vals for coords,vals in sj_dict[chrom].items()]):
                if args.FILTER:
                    motif = item[3]
                    if motif not in [1,2]:
                        continue
                
                print(chrom+'\t'+'\t'.join([str(i) for i in item]))
        elif args.FORMAT == 'bed':
            for item in sorted([tuple([int(i) for i in coords.split(':')])+vals for coords,vals in sj_dict[chrom].items()]):
                left, right, strand, motif, anno, unique, multi, overhang = item
                if args.FILTER:
                    if motif not in [1,2]:
                        continue
                
                if strand == 1:
                    s = '+'
                elif strand == 2:
                    s = '-'
                else:
                    s = '.'
                print('{}\t{}\t{}\t.\t{}\t{}'.format(
                    chrom,
                    left-1,
                    right,
                    unique,
                    s
                ))
    
