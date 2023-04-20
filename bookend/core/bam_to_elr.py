#!/usr/bin/python
# -*- coding: utf-8 -*-

from distutils import extension
import sys
import os
if __name__ == '__main__':
    sys.path.append('../../bookend')

from bookend.core.cython_utils._rnaseq_utils import RNAseqDataset
import pysam
from bookend.core.elr_sort import ELRsorter
from bookend.core.elr_to_bed import ELRtoBEDconverter

class BAMtoELRconverter:
    def __init__(self, args):
        """Parses input arguments for converting BAM to ELR"""
        print(args)
        self.source = args['SOURCE']
        self.genome = args['GENOME']
        self.reference = args['REFERENCE']
        self.splice = args['SPLICE']
        self.chrom_names = args['CHROM_NAMES']
        self.data_type = args['DATA_TYPE']
        self.stranded = args['STRANDED']
        self.reverse = args['REVERSE']
        self.stranded = self.stranded or self.reverse
        self.header = args['HEADER']
        self.start = args['START']
        self.capped = args['CAPPED']
        self.end = args['END']
        self.no_ends = args['NO_ENDS']
        self.secondary = args['SECONDARY']
        self.output = args['OUTPUT']
        self.start_seq = args['START_SEQ']
        self.end_seq = args['END_SEQ']
        self.untrimmed = args['UNTRIMMED']
        self.record_artifacts = args['RECORD_ARTIFACTS']
        self.mismatch_rate = args['MM_RATE']
        self.sj_shift = args['SJ_SHIFT']
        self.minlen_strict = args['MINLEN_STRICT']
        self.minlen_loose = args['MINLEN_LOOSE']
        self.input = args['INPUT']
        self.error_rate = args['ERROR_RATE']
        self.remove_noncanonical = args['REMOVE_NONCANONICAL']
        if self.start or self.end or self.capped:
            self.stranded = True
        
        self.ext = 'elr'
        if self.output is None:
            self.output = '{}.{}'.format(self.input,self.ext)
        elif '.' in self.output:
            self.ext = self.output.split('.')[-1]
        else:
            self.output = '{}.{}'.format(self.output,self.ext)
        
        if self.start_seq.lower == 'none':
            self.start_seq = ''
        
        if self.end_seq.lower == 'none':
            self.end_seq = ''
        
        if self.input.split('.')[-1].lower() not in ['bam','sam']:
            print("\nERROR: input file must be BAM/SAM format.")
            sys.exit(1)
        
        if self.output.split('.')[-1].lower() not in ['elr','bed', 'bed12']:
            print("\nERROR: output file must be ELR or BED format.")
            sys.exit(1)
        
        self.output_dict = {}
        self.tempout = '_unsorted.'+self.output
        self.tempout_file = open(self.tempout, 'w')
        if self.ext.lower() in ['bed','bed12']:
            self.output_format = 'bed'
        else:
            self.output_format = 'elr'

        self.config_dict = {
            'source':self.source,
            's_tag':self.start,
            'e_tag':self.end,
            'capped':self.capped,
            'stranded':self.stranded,
            'reverse':self.reverse,
            'start_seq':self.start_seq,
            'end_seq':self.end_seq,
            'minlen_strict':self.minlen_strict,
            'minlen_loose':self.minlen_loose,
            'mismatch_rate':self.mismatch_rate,
            'sj_shift':self.sj_shift,
            'remove_noncanonical':self.remove_noncanonical,
            'labels_are_trimmed':not self.untrimmed,
            'quality_filter':True,
            'reference':self.reference,
            'sj':self.splice,
        }
        if self.no_ends:
            self.config_dict['s_tag'] = False
            self.config_dict['e_tag'] = False
            self.config_dict['capped'] = False
            self.config_dict['start_seq'] = ''
            self.config_dict['end_seq'] = ''
        
        if self.data_type is not None:
            if self.data_type.upper() in ['ONT','NANOPORE','OXFORD']:
                """Reads are from Oxford Nanopore cDNA or PCR-cDNA kits, trimmed and oriented by bookend label or pychopper."""
                self.config_dict['stranded'] = True
                self.config_dict['max_headclip'] = 10
                if self.untrimmed:
                    self.config_dict['max_headclip'] = 120
                    self.config_dict['stranded'] = False
                    self.config_dict['s_tag'] = True
                    self.config_dict['e_tag'] = True
            elif self.data_type.upper() in ['PACBIO', 'ISOSEQ', 'ISOSEQ3', 'FLNC']:
                """Reads are PacBio FLNCs, downstream of lima."""
                self.config_dict['max_headclip'] = 4
                self.config_dict['s_tag'] = True
                self.config_dict['e_tag'] = True
                self.config_dict['quality_filter'] = True
            elif self.data_type.upper() in ['ONT-RNA','ONT_RNA','DIRECT_RNA', 'DIRECT-RNA']:
                """Reads are from Oxford Nanopore direct RNA kit, downstream of basecalling."""
                self.config_dict['stranded'] = True
                self.config_dict['labels_are_trimmed'] = False
                self.config_dict['quality_filter'] = True
            elif self.data_type.strip('0123456789').upper() in ['SMART','SMARTER','SMARTSEQ','SMART-SEQ']:
                """Reads are from a SMART protocol, labeled by bookend label."""
                self.config_dict['stranded'] = False
                self.config_dict['labels_are_trimmed'] = True
                self.config_dict['quality_filter'] = True
            else:
                print("\nERROR: --data_type not recognized.")
                print("Currently supported: ONT, PACBIO, DIRECT_RNA, SMARTSEQ")
                sys.exit(1)
        
        save = pysam.set_verbosity(0)
        self.bam_in = pysam.AlignmentFile(self.input)
        save = pysam.set_verbosity(save)
        if self.source is None:
            self.source = self.bam_in.header['PG'][0]['ID']
        
        self.dataset = RNAseqDataset(
            chrom_array=self.bam_in.header.references, 
            chrom_lengths=list(self.bam_in.header.lengths),
            source_array=[self.source],
            config=self.config_dict,
            genome_fasta=self.genome
        )
        if self.chrom_names is not None:
            # Swap chromosome names for the given 2-column TSV
            chromdict = {l.split('\t')[0]:l.rstrip().split('\t')[1] for l in open(self.chrom_names,'r')}
            chromdictrev = {v:k for k,v in chromdict.items()}
            for i in range(len(self.dataset.chrom_array)):
                c = self.dataset.chrom_array[i]
                if c in chromdict.keys():
                    self.dataset.genome[chromdict[c]] = self.dataset.genome[c]
                    self.dataset.chrom_array[i] = chromdict[c]
                elif c in chromdictrev.keys():
                    self.dataset.chrom_array[i] = chromdictrev[c]
                    self.dataset.genome[chromdictrev[c]] = self.dataset.genome[c]
        
        self.sort_args = {
            'OUT':self.output,
            'FORCE':True,
            'INPUT':self.tempout
        }
        if self.output_format == 'bed':
            self.sort_args['OUT'] += '.elr'
        
        self.failures = []
        self.generator = self.generate_bam_entries()
    
    def generate_bam_entries(self):
        """Group all BAM lines with the same read ID
        and yield them as a list of pysam objects"""
        current_ID = None
        bam_lines = []
        for line in self.bam_in:
            ID = line.query_name
            if ID == current_ID or current_ID is None:
                bam_lines.append(line)
                current_ID = ID
            else:
                yield bam_lines
                bam_lines = [line]
                current_ID = ID
        
        yield bam_lines

    def output_lines(self, lines, output):
        """Takes a list of bed lines and writes
        them to the output stream.
        """
        if output == 'stdout':
            for output_string in lines:
                print(output_string)
        else:
            output.write('\n'.join(lines)+'\n')

    def write_elr(self, output):
        out_strings = [mapping.write_as_elr(record_artifacts=self.record_artifacts).rstrip() for mapping in self.dataset.read_list]
        self.output_lines(out_strings, output)

    def process_entry(self, bam_lines):
        self.dataset.read_list = []
        # self.dataset.add_read_from_BAM(bam_lines, ignore_ends=self.no_ends, secondary=self.secondary, error_rate=self.error_rate)
        self.dataset.add_read_from_BAM(bam_lines)
        if len(self.dataset.read_list) > 0:
            self.write_elr(self.tempout_file)
        else:
            self.failures += bam_lines
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                          {}\n".format(self.input)
        options_string += "  Reference genome file:               {}\n".format(self.genome)
        options_string += "  Reference splice junction file:      {}\n".format(self.splice)
        options_string += "  Output file (-o):                    {}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Reads start at RNA 5' ends (-s):     {}\n".format(self.start)
        options_string += "  Reads are from capped RNA (-c):      {}\n".format(self.capped)
        options_string += "  Reads end at RNA 3' ends (-e):       {}\n".format(self.end)
        options_string += "  Untrimmed input reads (--untrimmed)  {}\n".format(self.untrimmed)
        options_string += "  Strand specific (--stranded):        {}\n".format(self.stranded)
        options_string += "  Strand specific reverse (--reverse): {}\n".format(self.reverse)
        options_string += "  Start tag suffix (--start_seq):      {}\n".format(self.start_seq)
        options_string += "  End tag prefix (--end_seq):          {}\n".format(self.end_seq)
        options_string += "  Adjust splice sites (--sj_shift):    {}\n".format(self.sj_shift)
        options_string += "  *** Filters ***\n"
        options_string += "  --record_artifacts:                  {}\n".format(self.record_artifacts)
        options_string += "  --mismatch_rate:                     {}\n".format(self.mismatch_rate)
        options_string += "  --remove_noncanonical:               {}\n".format(self.remove_noncanonical)
        options_string += "  Perfect minlen (--minlen_strict):    {}\n".format(self.minlen_strict)
        options_string += "  Relaxed minlen (--minlen_loose):     {}\n".format(self.minlen_loose)
        options_string += "  Secondary alignments (--secondary):  {}\n".format(self.secondary)

        if not self.genome:
            options_string += "\nWARNING: cap detection and artifact masking can only be done if a reference genome is provided."
            options_string += "\nProvide a genome fasta with --genome /path/to/fasta"
            if self.remove_noncanonical:
                options_string += "\nWARNING: noncanonical splice junctions can only be detected if --genome is provided."
        
        if not self.splice and not self.reference and self.sj_shift:
            options_string += "\nWARNING: splice junction correction can only be performed with a splice junction reference file."
            options_string += "\nProvide a BED12/GTF/GFF junction file with --reference /path/to/reference"
            options_string += "\nor a STAR/BED6 intron file with --splice /path/to/introns."
        
        return options_string
    
    def display_summary(self):
        if self.record_artifacts:
            # summary = 'len\tS\ts\tE\te\n'
            # max_len = max([
            #     max(self.dataset.label_tally['S']) if len(self.dataset.label_tally['S']) else 0,
            #     max(self.dataset.label_tally['s']) if len(self.dataset.label_tally['s']) else 0, 
            #     max(self.dataset.label_tally['E']) if len(self.dataset.label_tally['E']) else 0, 
            #     max(self.dataset.label_tally['e']) if len(self.dataset.label_tally['e']) else 0
            # ])
            # for i in range(max_len+1):
            #     summary += '{}\t{}\t{}\t{}\t{}\n'.format(
            #         i,
            #         self.dataset.label_tally['S'][i],
            #         self.dataset.label_tally['s'][i],
            #         self.dataset.label_tally['E'][i],
            #         self.dataset.label_tally['e'][i]
            #     )
            
            summary = 'Total\t{}\t{}\t{}\t{}\n'.format(
                sum(self.dataset.label_tally['S'].values()),
                sum(self.dataset.label_tally['s'].values()),
                sum(self.dataset.label_tally['E'].values()),
                sum(self.dataset.label_tally['e'].values())
            )
        else:
            # summary = 'len\tS\tE\n'
            # max_len = max([
            #     max(self.dataset.label_tally['S']) if len(self.dataset.label_tally['S']) else 0,
            #     max(self.dataset.label_tally['E']) if len(self.dataset.label_tally['E']) else 0
            # ])
            # for i in range(max_len+1):
            #     summary += '{}\t{}\t{}\n'.format(
            #         i,
            #         self.dataset.label_tally['S'][i],
            #         self.dataset.label_tally['E'][i]
            #     )
            
            summary = 'Total\t{}\t{}\n'.format(
                sum(self.dataset.label_tally['S'].values()),
                sum(self.dataset.label_tally['E'].values())
            )
            
        
        return summary
    
    def run(self):
        """Executes end labeling on all reads."""
        if self.output != 'stdout':
            print(self.display_options())
        
        self.output_lines(self.dataset.dump_header(),self.tempout_file)       
        for entry in self.generator:
            self.process_entry(entry)
        
        self.tempout_file.close()
        Sorter = ELRsorter(self.sort_args)
        Sorter.run()
        os.remove(self.tempout)
        if self.output_format == 'bed':
            convert_args = {
                'INPUT':self.sort_args['OUT'],
                'HEADER':None,
                'OUTPUT':self.output
            }
            Converter = ELRtoBEDconverter(convert_args)
            print('Converting to BED...')
            Converter.run()
            os.remove(self.sort_args['OUT'])
        
        print(self.display_summary())



if __name__ == '__main__':
    from argument_parsers import bam_to_elr_parser as parser
    args = vars(parser.parse_args())
    obj = BAMtoELRconverter(args)
    sys.exit(obj.run())

    # TESTING #
    # import cProfile
    # import pstats
    # profile = cProfile.Profile()
    # profile.runcall(obj.run)
    # ps = pstats.Stats(profile)
    # ps.print_stats()

