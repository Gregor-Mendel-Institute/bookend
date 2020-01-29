import re
import sys
import argparse
import _fasta_utils as fu
from collections import Counter
#import __end_label_old as el
if sys.version_info >= (3,0):
    izip = zip
else:
    from itertools import izip

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser(
    description="Trims and labels FASTQ reads that have paired 5' and 3' adapter sequences."
)
parser.add_argument(
    '-S', '--start', dest='START',
    help="Template switching primer, marks the RNA 5' terminus.",
    type=str, default='AAGCAGTGGTATCAACGCAGAGTACGGG'
)
parser.add_argument(
    '-E', '--end', dest='END',
    help="Reverse transcription primer, marks (reverse complement of) the RNA 3' terminus.",
    # type=str, default='AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTT+'
    type=str, default='TTTTTTTTTTTTTTTTTTTT+'
)
parser.add_argument(
    '--strand', dest='STRAND',
    help="Orientation of mate1 with respect to the RNA strand",
    type=str, default='unstranded', choices=['forward','reverse','unstranded']
)
parser.add_argument(
    '--out1', dest='OUT1',
    help="Destination file for trimmed and labeled mate 1 FASTQ file.",
    type=str, default='end_label.1.fastq'
)
parser.add_argument(
    '--out2', dest='OUT2',
    help="Destination file for trimmed and labeled mate 2 FASTQ file.",
    type=str, default='end_label.2.fastq'
)
parser.add_argument(
    '--solo', dest='SOLO_OUT',
    help="(paired-end only) Destination file for reads that are only informative single-end.",
    type=str, default='end_label.solo.fastq'
)
parser.add_argument(
    '--min_start', dest='MIN_START',
    help="Minimum number of nucleotides that must match to the start adapter sequence.",
    type=int, default=5
)
parser.add_argument(
    '--min_end', dest='MIN_END',
    help="Minimum number of nucleotides that must match to the end adapter sequence.",
    type=int, default=8
)
parser.add_argument(
    '--suppress_untrimmed', dest='SUPPRESS_UNTRIMMED',
    help="If no trimming occurred, do not write the read to output.",
    default=False, action='store_true'
)
parser.add_argument(
    '--verbose', dest='VERBOSE',
    help="Pause after each trimming decision and display result.",
    default=False, action='store_true'
)
parser.add_argument(
    '--minlen', dest='MINLEN',
    help="Minimum sequence length to keep.",
    type=int, default=18
)
parser.add_argument(
    '--mismatch_rate', dest='MM_RATE',
    help="Highest allow proportion of mismatches.",
    type=float, default=.06
)
parser.add_argument(
    '--minqual', dest='MINQUAL',
    help="Suppresses any trimmed sequences with lower than this mean phred quality score.",
    type=float, default=25
)
parser.add_argument(
    dest='FASTQ',
    help="Input FASTQ file(s). 1 for single-end, 2 for paired-end",
    type=str,
    nargs='+'
)
args = parser.parse_args()
if args.START.lower() == 'none':
    args.START = ''

if args.END.lower() == 'none':
    args.END = ''

#############
# FUNCTIONS #
#############

def display_trim(mate1, mate2, trim1, trim2, label):
    """Print the trimmed decision to stdout"""
    color = '2;30;47'
    if label.startswith('S'):
        color = '7;34;40'
    if label.startswith('E'):
        color = '6;30;41'
    if 'S' in label and 'E' in label:
        color = '4;30;45'
    print('{}\t{}'.format(mate1,mate2))
    print('\x1b[{}m{}\t{}\t{}\x1b[0m '.format(color,label,trim1,trim2))

def print_label_details(labeldict):
    """Given a Counter object of end labels, print a table of their frequency."""
    total = sum(labeldict.values())
    no_tag = labeldict['']
    s_tag = sum([labeldict[k] for k in labeldict.keys() if 'S' in k])
    e_tag = sum([labeldict[k] for k in labeldict.keys() if 'E' in k])
    print("Total reads processed: {}".format(total))
    print("Start tag: {} ({}%)".format(s_tag, int(s_tag/total*1000)/10))
    print("End tag:   {} ({}%)".format(e_tag, int(e_tag/total*1000)/10))
    print("Unlabeled: {} ({}%)".format(no_tag, int(no_tag/total*1000)/10))
    
    multilabel = sorted([k for k in labeldict.keys() if 'S' in k and 'E' in k])
    label_lengths = [int(i.strip('SE')) for i in labeldict.keys() if i != '' and i not in multilabel]
    if len(label_lengths) > 0:
        print("\nlen\tS_tag\tE_tag")
        max_len = max(label_lengths)
        for i in range(1, max_len+1):
            print('{}\t{}\t{}'.format(i, labeldict['S{}'.format(i)], labeldict['E{}'.format(i)]))
        
        if len(multilabel) > 0:
            print("Multi-label reads:")
            for m in multilabel:
                print('{}\t{}'.format(m, labeldict[m]))


###################################

if __name__ == '__main__':
    file1 = None
    file2 = None
    namesplit = re.compile('[ /\t]')
    minstart = args.MIN_START
    minend = args.MIN_END
    minlen = args.MINLEN
    minqual = args.MINQUAL
    strand = args.STRAND
    labeldict = Counter()
    if len(args.FASTQ) == 2:
        experiment_type = "PE"
        file1=open(args.FASTQ[0],'r')
        file2=open(args.FASTQ[1],'r')
    elif len(args.FASTQ) == 1:
        experiment_type = "SE"
        file1=open(args.FASTQ[0],'r')
    else:
        print("ERROR: please provide either 1 (single-end) or 2 (paired-end) FASTQ files.")
        sys.exit(1)

    S5string = args.START
    if len(S5string) > 1:
        if S5string[-1] == '+': # 3'-terminal monomer is specified
            S5monomer = fu.IUPACnum[S5string[-2]]
            S3monomer = fu.IUPACnum[fu.complement(S5string[-2])]
            S5string = S5string[:-1]
        else:
            S5monomer = S3monomer = -1
    else:
        S5monomer = S3monomer = -1
    
    S3string = fu.complement(S5string)

    E5string = args.END
    if len(E5string) > 1:
        if E5string[-1] == '+': # Monomer is specified
            E5monomer = fu.IUPACnum[E5string[-2]]
            E3monomer = fu.IUPACnum[fu.complement(E5string[-2])]
            E5string = E5string[:-1]
        else:
            E5monomer = E3monomer = -1
    else:
        E5monomer = E3monomer = -1
    
    E3string = fu.complement(E5string)
    
    S5array = fu.nuc_to_int(S5string, 'J'*len(S5string))
    S3array = fu.nuc_to_int(S3string, 'J'*len(S3string))
    E5array = fu.nuc_to_int(E5string, 'J'*len(E5string))
    E3array = fu.nuc_to_int(E3string, 'J'*len(E3string))
    
    linecounter=0
    outfile1=open(args.OUT1,'w')
    file1_read = None
    file2_read = None
    if file1 and file2:
        outfile_solo=open(args.SOLO_OUT,'w')
        outfile2=open(args.OUT2,'w')
        for line1, line2 in izip(file1, file2):
            linecounter+=1
            if linecounter % 4 == 1:
                if file1_read and file2_read:
                    # Execute paired-end trimming
                    trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
                        file1_read[1], file1_read[3],
                        file2_read[1], file2_read[3],
                        S5array, S5monomer,
                        S3array, S3monomer, 
                        E5array, E5monomer, 
                        E3array, E3monomer,
                        strand, minstart, minend, minlen, minqual
                    )
                    if args.VERBOSE:
                        display_trim(file1_read[1], file2_read[1], trim1, trim2, label)
                    
                    labeldict[label] += 1
                    if len(trim1) > 0:
                        if not args.SUPPRESS_UNTRIMMED or label != '':
                            if trim2 is None:
                                if len(trim1) >= args.MINLEN:
                                    # A solo read was chosen
                                    outfile_solo.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                                        file1_read[0],
                                        label,
                                        trim1,
                                        file1_read[2],
                                        qtrm1
                                    ))
                            else:
                                if len(trim1) >= args.MINLEN and len(trim2) >= args.MINLEN:
                                    # Write trimmed, oriented, paired-end reads
                                    outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                                        file1_read[0],
                                        label,
                                        trim1,
                                        file1_read[2],
                                        qtrm1
                                    ))
                                    outfile2.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                                        file2_read[0],
                                        label,
                                        trim2,
                                        file2_read[2],
                                        qtrm2
                                    ))
                
                name1 = re.split(namesplit,line1)[0].rstrip()
                name2 = re.split(namesplit,line2)[0].rstrip()
                if name1 != name2:
                    name1 = name1[:-2]
                    name2 = name2[:-2]
                
                assert name1 == name2, "ERROR: mate pair names do not match"
                file1_read = [name1]
                file2_read = [name2]
            else:
                file1_read.append(line1.rstrip())
                file2_read.append(line2.rstrip())
        
        # Execute paired-end trimming for the last read
        trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
            file1_read[1], file1_read[3],
            file2_read[1], file2_read[3],
            S5array, S5monomer,
            S3array, S3monomer, 
            E5array, E5monomer, 
            E3array, E3monomer,
            strand, minstart, minend, minlen, minqual
        )
        if args.VERBOSE:
            display_trim(file1_read[1], file2_read[1], trim1, trim2, label)
        
        labeldict[label] += 1
        if len(trim1) > 0:
            if not args.SUPPRESS_UNTRIMMED or label != '':
                if trim2 is None:
                    if len(trim1) >= args.MINLEN:
                        # A solo read was chosen
                        outfile_solo.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                            file1_read[0],
                            label,
                            trim1,
                            file1_read[2],
                            qtrm1
                        ))
                else:
                    if len(trim1) >= args.MINLEN and len(trim2) >= args.MINLEN:
                        # Write trimmed, oriented, paired-end reads
                        outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                            file1_read[0],
                            label,
                            trim1,
                            file1_read[2],
                            qtrm1
                        ))
                        outfile2.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                            file2_read[0],
                            label,
                            trim2,
                            file2_read[2],
                            qtrm2
                        ))
        
        file1.close()
        file2.close()
        outfile1.close()
        outfile2.close()
        outfile_solo.close()
        print_label_details(labeldict)
    elif file1:
        for line1 in file1:
            linecounter+=1
            if linecounter % 4 == 1:
                if file1_read:
                    # Execute single-end trimming
                    trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
                        file1_read[1], file1_read[3],
                        '', '',
                        S5array, S5monomer,
                        S3array, S3monomer, 
                        E5array, E5monomer, 
                        E3array, E3monomer,
                        strand, minstart, minend, minlen, minqual
                    )
                    if args.VERBOSE:
                        display_trim(file1_read[1], '', trim1, '', label)
                    
                    labeldict[label] += 1
                    if len(trim1) > 0:
                        if not args.SUPPRESS_UNTRIMMED or label != '':
                            if len(trim1) >= args.MINLEN:
                                # A solo read was chosen
                                outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                                    file1_read[0],
                                    label,
                                    trim1,
                                    file1_read[2],
                                    qtrm1
                                ))
                    
                name1 = re.split(namesplit,line1)[0].rstrip()
                file1_read = [name1]
            else:
                file1_read.append(line1.rstrip())
        
        # Execute single-end trimming for the last read
        trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
            file1_read[1], file1_read[3],
            '', '',
            S5array, S5monomer,
            S3array, S3monomer, 
            E5array, E5monomer, 
            E3array, E3monomer,
            strand, minstart, minend, minlen, minqual
        )
        if args.VERBOSE:
            display_trim(file1_read[1], '', trim1, '', label)
        
        labeldict[label] += 1
        if len(trim1) > 0:
            if not args.SUPPRESS_UNTRIMMED or label != '':
                if len(trim1) >= args.MINLEN:
                    # A solo read was chosen
                    outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(
                        file1_read[0],
                        label,
                        trim1,
                        file1_read[2],
                        qtrm1
                    ))
        
        file1.close()
        outfile1.close()
        print_label_details(labeldict)
    else:
        print("ERROR: no FASTQ files parsed.")
        sys.exit(1)
