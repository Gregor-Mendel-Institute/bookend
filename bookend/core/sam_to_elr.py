import sys
import argparse
from cython_utils import _rnaseq_utils as ru
import pysam

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    "--source", dest='SOURCE',
    help="Source of SAM lines.",
    default=None, type=str
)
parser.add_argument(
    "-F", "--fasta", dest='FASTA',
    help="Genome FASTA file",
    default=None, type=str
)
parser.add_argument(
    "--stranded", dest='STRANDED',
    help="The read is strand-specific",
    default=False, action="store_true"
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
    "--no_ends", dest='NO_ENDS', default=False, action='store_true',
    help="(overrides other options) Ignore all 5' and 3' end labels."
)
parser.add_argument(
    "--bed", dest='BED_OUT', default=False, action='store_true',
    help="Output a 15-column end labeled BED file."
)
parser.add_argument(
    "--secondary", dest='SECONDARY', default=False, action='store_true',
    help="Keep secondary alignments."
)
parser.add_argument(
    "--header", dest='HEADER', type=str, default=None,
    help="Filepath to write ELR header."
)
parser.add_argument(
    "-O", "--output", dest='OUTPUT', type=str, default='stdout',
    help="Filepath to write end-labeled file."
)
parser.add_argument(
    "--start_seq", dest='START_SEQ',
    help="Sequence of the oligo that marks a 5' read (sense)",
    default='ACGGG', type=str
)
parser.add_argument(
    "--end_seq", dest='END_SEQ',
    help="Sequence of the oligo that marks a 3' read (sense)",
    default='RRRRRRRRRRRRRRR', type=str
)
parser.add_argument(
    "--record_artifacts", dest='RECORD_ARTIFACTS', default=False, action='store_true',
    help="Reports artifact-masked S/E labels as lowercase s/e."
)
parser.add_argument(
    "--split", dest='SPLIT', default=False, action='store_true',
    help="Separate reads into different files by their multimapping number."
)
parser.add_argument(
    "--mismatch_rate", dest='MM',
    help="Mismatch tolerance of S label matches",
    default=0.20, type=float
)
parser.add_argument(
    "--sj_shift", dest='SJ_SHIFT',
    help="Shift up to this many bases to find a canonical splice junction",
    default=2, type=int
)
parser.add_argument(
    "--minlen", dest='MINLEN',
    help="Minimum number of nucleotides for an aligned read.",
    default=20, type=int
)
parser.add_argument(
    "INPUT", type=str, default=None, help="Input BAM/SAM file"
)
args = parser.parse_args()
if args.START or args.END or args.CAPPED:
    args.STRANDED = True

output_file = 'stdout'
output_dict = {}
if not args.SPLIT:
    if args.OUTPUT != 'stdout':
        output_file = open(args.OUTPUT, 'w')

if args.HEADER is None:
    header_file = output_file
elif args.HEADER != 'stdout':
    header_file = open(args.HEADER, 'w')
else:
    header_file = 'stdout'

if args.BED_OUT:
    output_format = 'bed'
else:
    output_format = 'elr'

config_dict = {
    'source':args.SOURCE,
    's_tag':args.START,
    'e_tag':args.END,
    'capped':args.CAPPED,
    'stranded':args.STRANDED,
    'start_seq':args.START_SEQ,
    'end_seq':args.END_SEQ,
    'minlen':args.MINLEN,
    'mismatch_rate':args.MM,
    'sj_shift':args.SJ_SHIFT
}
if args.NO_ENDS:
    config_dict['s_tag'] = False
    config_dict['e_tag'] = False
    config_dict['capped'] = False
    config_dict['start_seq'] = ''
    config_dict['end_seq'] = ''

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


def write_read(mapping_list, chrom_array, source_array, output, minlen=args.MINLEN, bed=args.BED_OUT, record_artifacts=args.RECORD_ARTIFACTS):
    nmap = len(mapping_list)
    for i,mapping in enumerate(mapping_list):
        match_length = mapping.get_length()
        if match_length >= minlen:
            if bed:
                out_fields = mapping.write_as_bed(chrom_array, source_array, as_string=False, record_artifacts=record_artifacts)
                out_fields[6] = nmap
                out_fields[7] = i+1
                out_string = '\t'.join([str(f) for f in out_fields])
            else:
                out_string = mapping.write_as_elr(record_artifacts=record_artifacts)
            
            output_lines([out_string], output)

######################
# PROCESS INPUT FILE # 
######################
if __name__ == '__main__':
    if args.INPUT:
        if args.INPUT.split('.')[-1].lower() not in ['bam','sam']:
            print("\nERROR: input file must be BAM/SAM format.")
            parser.print_help()
            sys.exit(1)
        
        bam_in = pysam.AlignmentFile(args.INPUT)
    else:
        print("\nERROR: requires BAM/SAM file as input.")
        parser.print_help()
        sys.exit(1)
    
    source = ''
    if args.SOURCE is None:
        source = bam_in.header['PG'][0]['ID']
    else:
        source = args.SOURCE
    
    if not args.FASTA:
        dataset = ru.RNAseqDataset(
                chrom_array=bam_in.header.references, 
                chrom_lengths=list(bam_in.header.lengths), 
                source_array=[source],
                config=config_dict)
    else:
        dataset = ru.RNAseqDataset(
                source_array=[source], 
                genome_fasta=args.FASTA,
                config=config_dict)
    
    output_lines(dataset.dump_header(), header_file)
    
    current_ID    = None
    current_read  = [None]
    current_lines = []
    for line in bam_in:
        ID = line.query_name
        if ID == current_ID or current_ID is None:
            current_lines.append(line)
            current_ID = ID
        else:
            # A new read ID was encountered. Perform actions on current_read
            # and begin a new current_read with the new ID.
            #if current_ID == 'HWI-ST1253F_0152:2:1101:1229:73510#19608_TAAGGCGA_TAG=E11':
            #    sys.exit(1)
                
            dataset.add_read_from_BAM(current_lines, ignore_ends=args.NO_ENDS, secondary=args.SECONDARY)
            if args.SPLIT: # Separate reads by their number of mappings to the genome
                mm_num = len(dataset.read_list)
                if mm_num not in output_dict.keys():
                    # Start a new file when a new multimapping number is found
                    output_dict[mm_num] = open('{}.mm{}.elr'.format(args.SOURCE, mm_num), 'w')
                    # Add the header to all output files
                    output_lines(dataset.dump_header(), output_dict[mm_num])

                output_file = output_dict[mm_num]
            
            write_read(dataset.read_list, dataset.chrom_array, dataset.source_array, output_file)
            # sys.exit(0)
            dataset.read_list = []
            current_lines = [line]
            current_ID = ID
    
    dataset.add_read_from_BAM(current_lines, ignore_ends=args.NO_ENDS, secondary=args.SECONDARY)
    if args.SPLIT: # Separate reads by their number of mappings to the genome
        mm_num = len(dataset.read_list)
        if mm_num not in output_dict.keys():
            output_dict[mm_num] = open('{}.mm{}.elr'.format(args.SOURCE, mm_num), 'w')
            output_lines(dataset.dump_header(), output_dict[mm_num])
        
        output_file = output_dict[mm_num]
    
    write_read(dataset.read_list, dataset.chrom_array, dataset.source_array, output_file)
    if args.OUTPUT != 'stdout':
        output_file.close()

    if args.RECORD_ARTIFACTS:
        max_len = max([
            max(dataset.label_tally['S']),
            max(dataset.label_tally['s']), 
            max(dataset.label_tally['E']), 
            max(dataset.label_tally['e'])
        ])
        print('len\tS\ts\tE\te')
        for i in range(max_len+1):
            print('{}\t{}\t{}\t{}\t{}'.format(
                i,
                dataset.label_tally['S'][i],
                dataset.label_tally['s'][i],
                dataset.label_tally['E'][i],
                dataset.label_tally['e'][i]
            ))
