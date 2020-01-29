'''
Merges any number of 'SJ.out.tab' files (output from the STAR aligner),
or any number of 'SJ.bed' files (--sj_out output from sam_calculate_coverages.py).
Outputs to stdout.
Description of SJ.out.tab, from Star Manual 2.5.1a, copyright Alexander Dobin 2016:
    SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that
    STAR defines the junction start/end as intronic bases, while many other software define them as
    exonic bases. The columns have the following meaning:
        column 1: chromosome
        column 2: first base of the intron (1-based)
        column 3: last base of the intron (1-based)
        column 4: strand (0: undefined, 1: +, 2: -)
        column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
        AT/AC, 6: GT/AT
        column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
        column 7: number of uniquely mapping reads crossing the junction
        column 8: number of multi-mapping reads crossing the junction
        column 9: maximum spliced alignment overhang

Description of SJ.bed:
    column 1: chromosome
    column 2: leftmost base of intron
    column 3: rightmost base of intron
    column 4: name (not used)
    column 5: hit-normalized read depth
    column 6: strand (+, -, or .)

File type is inferred from the number of columns
in the first line of the first file. Do not mix file types.
'''
import sys
import argparse
from argparse import RawTextHelpFormatter

desc = (
    "Combines multiple SJ.out.tab or SJ.bed files."
)
# initilize argumentparser to read in commands
parser = argparse.ArgumentParser(
    description=desc,
    formatter_class=RawTextHelpFormatter
)
parser.add_argument("--min_unique", dest='MIN_UNIQUE', 
    help="Filter SJs with fewer unique reads.",
    default=0, type=int
)
parser.add_argument("--min_reps", dest='MIN_REPS', 
    help="Filter SJs detected in fewer than this many files.",
    default=1, type=int
)
parser.add_argument("--new", dest='NEW', 
    help="Keep only SJs not present in the reference SJDB",
    default=False, action='store_true'
)
parser.add_argument(
    "FILES", nargs='+'
)
args = parser.parse_args()

sj_dict = {}
filetype = None
for sj_file in args.FILES:
    for line in open(sj_file):
        l = line.rstrip().split('\t')
        if not filetype:
            if len(l) == 9:
                filetype = 'star'
            elif len(l) == 6:
                filetype = 'bed'
            else:
                print("ERROR: couldn't determine file type")
                sys.exit(1)
        if filetype == 'star':
            if args.NEW:
                if l[5]=='1':
                    continue

            chrom = l[0]
            coords = ':'.join(l[1:6])
            unique,multiple,overhang = l[6:]
            if chrom not in sj_dict:
                sj_dict[chrom] = {}
            x = sj_dict[chrom].get(coords,(0,0,0,0))
            y = (int(unique),int(multiple),int(overhang))
            sj_dict[chrom][coords] = (x[0]+y[0],x[1]+y[1],max([x[2],y[2]]),x[3]+1)
        elif filetype == 'bed':
            chrom = l[0]
            if l[5] == '+':
                strand = '1'
            elif l[5] == '-':
                strand = '2'
            else:
                strand = '0'
            
            coords = ':'.join([l[1],l[2],strand])
            score = float(l[4])
            if chrom not in sj_dict:
                sj_dict[chrom] = {}
            sj_dict[chrom][coords] = sj_dict[chrom].get(coords,float(0)) + score
            

for chrom in sorted(list(sj_dict.keys())):
    if filetype == 'star':
        for item in sorted([tuple([int(i) for i in coords.split(':')])+vals for coords,vals in sj_dict[chrom].items()]):
            count = item[-1]
            unique = item[-4]
            if count >= args.MIN_REPS and unique >= args.MIN_UNIQUE:
                print(chrom+'\t'+'\t'.join([str(i) for i in item[:-1]]))
    elif filetype == 'bed':
        for left,right,strand,score in sorted([tuple([int(i) for i in coords.split(':')]+[vals]) for coords,vals in sj_dict[chrom].items()]):
            if strand == 1:
                s = '+'
            elif strand == 2:
                s = '-'
            else:
                s = '.'
            print('{}\t{}\t{}\t.\t{}\t{}'.format(
                chrom,
                left,
                right,
                score,
                s
            ))
