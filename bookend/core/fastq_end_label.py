#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import sys
import gzip
if __name__ == '__main__':
    sys.path.append('../../bookend')

from collections import Counter
import bookend.core.cython_utils._fasta_utils as fu
if sys.version_info >= (3,0):
    izip = zip
else:
    from itertools import izip

class EndLabeler:
    def __init__(self, args):
        """Parses input arguments for end labeling FASTQ file(s)"""
        self.namesplit = re.compile('[ /\t]')
        self.input = args['FASTQ']
        self.file1 = None
        self.file2 = None
        self.verbose = args['VERBOSE']
        self.minstart = args['MIN_START']
        self.minend = args['MIN_END']
        self.minlen = args['MINLEN']
        self.minqual = args['MINQUAL']
        self.qualmask = args['QUALMASK']
        self.strand = args['STRAND']
        self.discard_untrimmed = args['DISCARD_UNTRIMMED']
        self.out1 = args['OUT1']
        self.out2 = args['OUT2']
        self.pseudomates = args['PSEUDOMATES']
        self.umi = args['UMI']
        self.umi_range = (0,0)
        if self.pseudomates:
            self.single_out = None
        else:
            self.single_out = args['SINGLE_OUT']
        
        self.mm_rate = args['MM_RATE']
        self.labeldict = Counter()
        self.open_input_files()
       
        self.s_label = '' if args['START'].lower() == 'none' else args['START']
        self.S5string = self.s_label
        if len(self.S5string) > 1:
            if self.umi == 'S':
                self.umi_range = (self.S5string.find('N')-len(self.S5string), self.S5string.rfind('N')+1-len(self.S5string))
                if self.umi_range[1]-self.umi_range[0] == 0:
                    print("ERROR: No string of N's found for UMI in S label.")
                    sys.exit(1)
            
            if self.S5string[-1] == '+': # 3'-terminal monomer is specified
                self.S5monomer = fu.IUPACnum[self.S5string[-2]]
                self.S3monomer = fu.IUPACnum[fu.complement(self.S5string[-2])]
                self.S5string = self.S5string[:-1]
            else:
                self.S5monomer = self.S3monomer = -1
        else:
            self.S5monomer = self.S3monomer = -1
        
        self.S3string = fu.complement(self.S5string)
        self.e_label = '' if args['END'].lower() == 'none' else args['END']
        self.E5string = self.e_label
        if len(self.E5string) > 1:
            if self.umi == 'E':
                self.umi_range = (self.E5string.find('N')-len(self.E5string), self.E5string.rfind('N')+1-len(self.S5string))
                if self.umi_range[1]-self.umi_range[0] == 0:
                    print("ERROR: No string of N's found for UMI in E label.")
                    sys.exit(1)
            
            if self.E5string[-1] == '+': # Monomer is specified
                self.E5monomer = fu.IUPACnum[self.E5string[-2]]
                self.E3monomer = fu.IUPACnum[fu.complement(self.E5string[-2])]
                self.E5string = self.E5string[:-1]
            else:
                self.E5monomer = self.E3monomer = -1
        else:
            self.E5monomer = self.E3monomer = -1
        
        self.E3string = fu.complement(self.E5string)
        self.S5array = fu.nuc_to_int(self.S5string, 'J'*len(self.S5string))
        self.S3array = fu.nuc_to_int(self.S3string, 'J'*len(self.S3string))
        self.E5array = fu.nuc_to_int(self.E5string, 'J'*len(self.E5string))
        self.E3array = fu.nuc_to_int(self.E3string, 'J'*len(self.E3string))
        if self.single_out is not None:
            self.outfile_single=open(self.single_out,'w')
        
        if self.file2 is None: # Single-end input
            self.out1 = self.outfile1 = self.out2 = self.outfile2 =  None
        else: # Paired-end input
            self.outfile1 = open(self.out1,'w')
            self.outfile2 = open(self.out2,'w')
        
        self.fastq_generator = self.generate_fastq_entries()

    def run(self):
        """Executes end labeling on all reads."""
        print(self.display_options())
        for fastq_entry in self.fastq_generator:
            self.label_fastq_entry(fastq_entry)
        
        if self.single_out is not None:
            self.outfile_single.close()
        
        if self.file2 is not None:
            self.outfile1.close()
            self.outfile2.close()
        
        print(self.display_label_summary())
    
    def open_input_files(self):
        self.experiment_type = "SE"
        if self.input[0].endswith('gz'):
            self.file1 = gzip.open(self.input[0], 'rt')
        else:
            self.file1 = open(self.input[0],'r')
        
        if len(self.input) == 2:
            self.experiment_type = "PE"
            if self.input[1].endswith('gz'):
                self.file2 = gzip.open(self.input[1], 'rt')
            else:
                self.file2 = open(self.input[1],'r')
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend label |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file(s):                     {}\n".format(self.input)
        options_string += "  Output single-end (--single_out):  {}\n".format(self.single_out)
        options_string += "  Output paired-end mate 1 (--out1): {}\n".format(self.out1)
        options_string += "  Output paired-end mate 2 (--out2): {}\n".format(self.out2)
        options_string += "  *** Experiment parameters ***\n"
        options_string += "  Start tag (-S/--start):            {}\n".format(self.s_label)
        options_string += "  End tag (-E/--end):                {}\n".format(self.e_label)
        options_string += "  cDNA strand (--strand):            {}\n".format(self.strand)
        options_string += "  *** Filters ***\n"
        options_string += "  --discard_untrimmed:              {}\n".format(self.discard_untrimmed)
        options_string += "  --mismatch_rate:                   {}\n".format(self.mm_rate)
        options_string += "  S tag length min (--min_start):    {}\n".format(self.minstart)
        options_string += "  E tag length min (--min_end):      {}\n".format(self.minend)
        options_string += "  Trimmed length min (--minlen):     {}\n".format(self.minlen)
        options_string += "  Min avg. phred score (--minqual):  {}\n".format(self.minqual)
        options_string += "  Phred to mask as N (--qualmask):   {}\n".format(self.qualmask)
        options_string += "  Write pseudomates (--pseudomates): {}\n".format(self.pseudomates)
        return options_string

    def display_trim(self, mate1, mate2, trim1, trim2, label):
        """Returns a summary string of the trimming operation"""
        if mate2:
            trimstring = '[{}]\t{}|{} --> {}|{}'.format(label, mate1, mate2, trim1, trim2)
        else:
            trimstring = '[{}]\t{} --> {}'.format(label, mate1, trim1)
        
        return trimstring
    
    def display_label_summary(self):
        """Given a Counter object of end labels, returns a summary of their length frequency distribution."""
        total = sum(self.labeldict.values())
        if total == 0:
            summary_string = "No reads processed.\n"
        else:
            no_tag = self.labeldict['']+self.labeldict['XD']
            s_tag = sum([self.labeldict[k] for k in self.labeldict.keys() if 'S' in k])
            e_tag = sum([self.labeldict[k] for k in self.labeldict.keys() if 'E' in k])
            xl_tag = self.labeldict.get('XL', 0)
            xq_tag = self.labeldict.get('XQ', 0)
            kept = total-(xl_tag+xq_tag)
            summary_string = "Total reads processed: {}\n".format(total)
            summary_string += "Start labeled:   {} ({}%)\n".format(s_tag, int(s_tag/total*1000)/10)
            summary_string += "End labeled:     {} ({}%)\n".format(e_tag, int(e_tag/total*1000)/10)
            summary_string += "Unlabeled:       {} ({}%)\n".format(no_tag, int(no_tag/total*1000)/10)
            summary_string += "Removed (short): {} ({}%)\n".format(xl_tag, int(xl_tag/total*1000)/10)
            summary_string += "Removed (qual):  {} ({}%)\n".format(xq_tag, int(xq_tag/total*1000)/10)
            summary_string += "Total output:    {} ({}%)\n".format(kept, int(kept/total*1000)/10)

            multilabel = set([k for k in self.labeldict.keys() if 'S' in k and 'E' in k])
            for m in multilabel:
                s, e = m.lstrip('S').split('E')
                self.labeldict['S{}'.format(s)] += self.labeldict[m]
                self.labeldict['E{}'.format(e)] += self.labeldict[m]
            
            label_lengths = [int(i.strip('SE')) for i in self.labeldict.keys() if i != '' and i not in multilabel and 'X' not in i]
            if len(label_lengths) > 0:
                summary_string += "\nlen\tS_tag\tE_tag\n"
                max_len = max(label_lengths)
                for i in range(1, max_len+1):
                    s_count = self.labeldict['S{}'.format(i)]
                    e_count = self.labeldict['E{}'.format(i)]
                    if s_count > 0 or e_count > 0:
                        summary_string += '{}\t{}\t{}\n'.format(i, s_count, e_count)
                
                if len(multilabel) > 0:
                    summary_string += "Multi-label reads: {}".format(sum([self.labeldict[m] for m in multilabel]))
        
        return summary_string
    
    def generate_fastq_entries(self):
        """Generator object to iterate over one or two FASTQ files
        and package each read as a tuple"""
        file1_read = file2_read = None
        linecounter=0
        if self.file1 and self.file2:
            for line1, line2 in izip(self.file1, self.file2):
                linecounter+=1
                if linecounter % 4 == 1:
                    if file1_read and file2_read:
                        yield (file1_read, file2_read)
                    
                    name1 = re.split(self.namesplit,line1)[0].rstrip()
                    name2 = re.split(self.namesplit,line2)[0].rstrip()
                    if name1 != name2:
                        name1 = name1[:-2]
                        name2 = name2[:-2]
                    
                    assert name1 == name2, "ERROR: mate pair names do not match"
                    file1_read = [name1]
                    file2_read = [name2]
                else:
                    file1_read.append(line1.rstrip())
                    file2_read.append(line2.rstrip())
            
            yield (file1_read, file2_read)
        elif self.file1:
            for line1 in self.file1:
                linecounter+=1
                if linecounter % 4 == 1:
                    if file1_read:
                        yield (file1_read, None)
                    
                    name1 = re.split(self.namesplit,line1)[0].rstrip()
                    file1_read = [name1]
                else:
                    file1_read.append(line1.rstrip())
            
            yield (file1_read, None)
    
    def label_fastq_entry(self, fastq_entry):
        file1_read, file2_read = fastq_entry
        if file1_read and file2_read:
            # Execute paired-end trimming
            trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
                file1_read[1], file1_read[3], file2_read[1], file2_read[3],
                self.S5array, self.S5monomer, self.S3array, self.S3monomer, 
                self.E5array, self.E5monomer, self.E3array, self.E3monomer,
                self.strand, self.minstart, self.minend, self.minlen, self.minqual, self.qualmask, self.mm_rate, self.umi, self.umi_range
            )
            if self.verbose:
                print(self.display_trim(file1_read[1], file2_read[1], trim1, trim2, label))
            
            if len(trim1) > 0:
                if not self.discard_untrimmed or label != '':
                    if trim2 is None:
                        if len(trim1) >= self.minlen: # A single read was written from the pair
                            if fu.is_homopolymer(trim1):
                                self.labeldict['XQ'] += 1
                            else:
                                if self.pseudomates:
                                    self.outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file1_read[0], label, trim1, file1_read[2], qtrm1))
                                    self.outfile2.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file1_read[0], label, fu.rc(trim1), file1_read[2], qtrm1[::-1]))
                                else:
                                    self.outfile_single.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file1_read[0], label, trim1, file1_read[2], qtrm1))
                                
                                self.labeldict[label.split('_UMI=')[0]] += 1
                        else:
                            self.labeldict['XL'] += 1
                    else:
                        if len(trim1) >= self.minlen and len(trim2) >= self.minlen: # Write trimmed, oriented, paired-end reads
                            if fu.is_homopolymer(trim1) or fu.is_homopolymer(trim2):
                                self.labeldict['XQ'] += 1
                            else:
                                self.outfile1.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file1_read[0], label, trim1, file1_read[2], qtrm1))
                                self.outfile2.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file2_read[0], label, trim2, file2_read[2], qtrm2))
                                self.labeldict[label.split('_UMI=')[0]] += 1
                        else:
                            self.labeldict['XL'] += 1
                else:
                    self.labeldict['XD'] += 1
            else:
                if label == '':
                    self.labeldict['XQ'] += 1
                else:
                    self.labeldict['XL'] += 1
        elif file1_read: # Execute single-end trimming
            trim1, qtrm1, trim2, qtrm2, label = fu.terminal_trim(
                file1_read[1], file1_read[3], '', '',
                self.S5array, self.S5monomer, self.S3array, self.S3monomer, 
                self.E5array, self.E5monomer, self.E3array, self.E3monomer,
                self.strand, self.minstart, self.minend, self.minlen, self.minqual, self.qualmask, self.mm_rate, self.umi, self.umi_range
            )
            if self.verbose:
                print(self.display_trim(file1_read[1], '', trim1, '', label))
            
            if len(trim1) > 0:
                if fu.is_homopolymer(trim1):
                    self.labeldict['XQ'] += 1
                elif not self.discard_untrimmed or label != '':
                    if len(trim1) >= self.minlen: # A single read was output
                        self.outfile_single.write('{}_TAG={}\n{}\n{}\n{}\n'.format(file1_read[0], label, trim1, file1_read[2], qtrm1))
                        self.labeldict[label.split('_UMI=')[0]] += 1
                    else:
                        self.labeldict['XL'] += 1
                else:
                    self.labeldict['XD'] += 1
            else:
                if label == '':
                    self.labeldict['XQ'] += 1
                else:
                    self.labeldict['XL'] += 1

##########################
if __name__ == '__main__':
    from argument_parsers import end_label_parser
    args = vars(end_label_parser.parse_args())  
    obj = EndLabeler(args)
    sys.exit(obj.run())
