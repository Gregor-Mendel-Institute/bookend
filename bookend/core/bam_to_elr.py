#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
if __name__ == '__main__':
    sys.path.append('../../bookend')

import bookend.core.cython_utils._rnaseq_utils as ru
import pysam

class BAMtoELRconverter:
    def __init__(self, args):
        """Parses input arguments for converting BAM to ELR"""
        self.source = args['SOURCE']
        self.genome = args['GENOME']
        self.stranded = args['STRANDED']
        self.header = args['HEADER']
        self.start = args['START']
        self.capped = args['CAPPED']
        self.end = args['END']
        self.no_ends = args['NO_ENDS']
        self.bed_out = args['BED_OUT']
        self.secondary = args['SECONDARY']
        self.output = args['OUTPUT']
        self.start_seq = args['START_SEQ']
        self.end_seq = args['END_SEQ']
        self.record_artifacts = args['RECORD_ARTIFACTS']
        self.split = args['SPLIT']
        self.mismatch_rate = args['MM_RATE']
        self.sj_shift = args['SJ_SHIFT']
        self.minlen_strict = args['MINLEN_STRICT']
        self.minlen_loose = args['MINLEN_LOOSE']
        self.input = args['INPUT']
        self.error_rate = args['ERROR_RATE']
        self.remove_noncanonical = args['REMOVE_NONCANONICAL']
        if self.start or self.end or self.capped:
            self.stranded = True
        
        if self.output is None:
            self.output = self.input+'.elr'
        
        if self.start_seq.lower == 'none':
            self.start_seq = ''
        
        if self.end_seq.lower == 'none':
            self.end_seq = ''
        
        self.output_file = 'stdout'
        self.output_dict = {}
        if not self.split:
            if self.output != 'stdout':
                self.output_file = open(self.output, 'w')

        if self.header is None:
            self.header_file = self.output_file
        elif self.header != 'stdout':
            self.header_file = open(self.header, 'w')
        else:
            self.header_file = 'stdout'

        if self.bed_out:
            self.output_format = 'bed'
        else:
            self.output_format = 'elr'

        self.config_dict = {
            'source':self.source,
            's_tag':self.start,
            'e_tag':self.end,
            'capped':self.capped,
            'stranded':self.stranded,
            'start_seq':self.start_seq,
            'end_seq':self.end_seq,
            'minlen_strict':self.minlen_strict,
            'minlen_loose':self.minlen_loose,
            'mismatch_rate':self.mismatch_rate,
            'sj_shift':self.sj_shift,
            'remove_noncanonical':self.remove_noncanonical
        }
        if self.no_ends:
            self.config_dict['s_tag'] = False
            self.config_dict['e_tag'] = False
            self.config_dict['capped'] = False
            self.config_dict['start_seq'] = ''
            self.config_dict['end_seq'] = ''
        
        if self.input.split('.')[-1].lower() not in ['bam','sam']:
            print("\nERROR: input file must be BAM/SAM format.")
            sys.exit(1)
        
        self.bam_in = pysam.AlignmentFile(self.input)
        if self.source is None:
            self.source = self.bam_in.header['PG'][0]['ID']
        
        self.dataset = ru.RNAseqDataset(
            chrom_array=self.bam_in.header.references, 
            chrom_lengths=list(self.bam_in.header.lengths),
            source_array=[self.source],
            config=self.config_dict,
            genome_fasta=self.genome
        )
        
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

    def write_reads(self, output):
        out_strings = [mapping.write_as_elr(record_artifacts=self.record_artifacts).rstrip() for mapping in self.dataset.read_list]
        self.output_lines(out_strings, output)

    def process_entry(self, bam_lines):
        self.dataset.read_list = []
        self.dataset.add_read_from_BAM(bam_lines, ignore_ends=self.no_ends, secondary=self.secondary, error_rate=self.error_rate)
        if self.split: # Separate reads by their number of mappings to the genome
            mm_num = len(self.dataset.read_list)
            if mm_num not in self.output_dict.keys():
                # Start a new file when a new multimapping number is found
                self.output_dict[mm_num] = open('{}.mm{}.elr'.format(self.source, mm_num), 'w')
                # Add the header to all output files
                self.output_lines(self.dataset.dump_header(), self.output_dict[mm_num])

            output_file = self.output_dict[mm_num]
        else:
            output_file = self.output_file
        
        if len(self.dataset.read_list) > 0:
            self.write_reads(output_file)
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend make-elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:                         {}\n".format(self.input)
        options_string += "  Reference genome file:              {}\n".format(self.genome)
        options_string += "  Output file (-o):                   {}\n".format(self.output)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Reads start at RNA 5' ends (-s):    {}\n".format(self.start)
        options_string += "  Reads are from capped RNA (-c):     {}\n".format(self.capped)
        options_string += "  Reads end at RNA 3' ends (-e):      {}\n".format(self.end)
        options_string += "  Strand specific (--stranded):       {}\n".format(self.stranded)
        options_string += "  Start tag suffix (--start_seq):     {}\n".format(self.start_seq)
        options_string += "  End tag prefix (--end_seq):         {}\n".format(self.end_seq)
        options_string += "  Adjust splice sites (--sj_shift):   {}\n".format(self.sj_shift)
        options_string += "  *** Filters ***\n"
        options_string += "  --record_artifacts:                 {}\n".format(self.record_artifacts)
        options_string += "  --mismatch_rate:                    {}\n".format(self.mismatch_rate)
        options_string += "  Perfect alignment minlen (--minlen_strict): {}\n".format(self.minlen_strict)
        options_string += "  Relaxed alignment minlen (--minlen_loose):  {}\n".format(self.minlen_loose)
        options_string += "  Secondary alignments (--secondary): {}\n".format(self.secondary)
        return options_string
    
    def display_summary(self):
        if self.record_artifacts:
            summary = 'len\tS\ts\tE\te\n'
            max_len = max([
                max(self.dataset.label_tally['S']) if len(self.dataset.label_tally['S']) else 0,
                max(self.dataset.label_tally['s']) if len(self.dataset.label_tally['s']) else 0, 
                max(self.dataset.label_tally['E']) if len(self.dataset.label_tally['E']) else 0, 
                max(self.dataset.label_tally['e']) if len(self.dataset.label_tally['e']) else 0
            ])
            for i in range(max_len+1):
                summary += '{}\t{}\t{}\t{}\t{}\n'.format(
                    i,
                    self.dataset.label_tally['S'][i],
                    self.dataset.label_tally['s'][i],
                    self.dataset.label_tally['E'][i],
                    self.dataset.label_tally['e'][i]
                )
            
            summary += 'Total\t{}\t{}\t{}\t{}\n'.format(
                sum(self.dataset.label_tally['S'].values()),
                sum(self.dataset.label_tally['s'].values()),
                sum(self.dataset.label_tally['E'].values()),
                sum(self.dataset.label_tally['e'].values())
            )
        else:
            summary = 'len\tS\tE\n'
            max_len = max([
                max(self.dataset.label_tally['S']) if len(self.dataset.label_tally['S']) else 0,
                max(self.dataset.label_tally['E']) if len(self.dataset.label_tally['E']) else 0
            ])
            for i in range(max_len+1):
                summary += '{}\t{}\t{}\n'.format(
                    i,
                    self.dataset.label_tally['S'][i],
                    self.dataset.label_tally['E'][i]
                )
            
            summary += 'Total\t{}\t{}\n'.format(
                sum(self.dataset.label_tally['S'].values()),
                sum(self.dataset.label_tally['E'].values())
            )
            
        
        return summary
    
    def run(self):
        """Executes end labeling on all reads."""
        print(self.display_options())
        self.output_lines(self.dataset.dump_header(),self.output_file)
        for entry in self.generator:
            self.process_entry(entry)
        
        if self.output_file != 'stdout':
            self.output_file.close()
        
        for mm_out in self.output_dict.values():
            mm_out.close()
        
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

