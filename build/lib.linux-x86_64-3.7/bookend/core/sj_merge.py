#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
if __name__ == '__main__':
    sys.path.append('../../bookend')

class SJobject:
    def __init__(self, line, linetype=None):
        """Parses a line of SJ.out.tab or SJ.bed to an SJ object.
        O-indexed open coordinates."""
        fields = line.rstrip().split('\t')
        if linetype is None:
            if len(fields) == 9:
                linetype = 'star'
            else:
                linetype = 'bed'
        
        self.count = 1
        self.linetype = linetype
        self.chrom = fields[0]
        self.end = int(fields[2])
        if linetype == 'star':
            self.start = int(fields[1]) - 1
            self.strand = ['.', '+', '-'][int(fields[3])]
            self.motif = fields[4]
            self.new = fields[5] == '0'
            self.unique = int(fields[6])
            self.multi = int(fields[7])
            self.overhang = int(fields[8])
            self.name = None
        else:
            self.start = int(fields[1])
            self.strand = fields[5]
            self.motif = 0
            self.new = True
            self.unique = int(fields[4])
            self.multi = 0
            self.overhang = 0
            self.name = fields[3]
        
        self.hash = '{}_{}_{}_{}'.format(self.chrom, self.start, self.end, self.strand)
        self.comparator = (self.chrom, self.start, self.end)
    
    def __eq__(self, other): return self.comparator == other.comparator
    def __gt__(self, other): return self.comparator > other.comparator
    def __ge__(self, other): return self.comparator >= other.comparator
    def __lt__(self, other): return self.comparator < other.comparator
    def __le__(self, other): return self.comparator <= other.comparator
    def __ne__(self, other): return self.comparator != other.comparator
    
    def merge(self, other):
        """Combines the information of two compatible SJ entries"""
        if self.hash == other.hash:
            self.count += other.count
            self.unique += other.unique
            self.multi += other.multi
            self.overhang = max(self.overhang, other.overhang)
            if self.name is None:
                self.name = other.name
    
    def as_string(self, linetype='star'):
        """Writes the SJobject as a SJ.out.tab or SJ.bed format"""
        if linetype == 'star':
            return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                self.chrom, self.start+1, self.end, 
                {'.':0, '+':1, '-':2}[self.strand],
                self.motif,
                0 if self.new else 1,
                self.unique,
                self.multi,
                self.overhang
            )
        else:
            return '{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                self.chrom, self.start, self.end, self.name, self.unique+self.multi, self.strand
            )

class SJmerger:
    def __init__(self, args):
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
        in the first line of the first file.
        '''
        self.input = args['INPUT']
        self.output = args['OUT']
        self.min_unique = args['MIN_UNIQUE']
        self.min_reps = args['MIN_REPS']
        self.format = args['FORMAT']
        self.new = args['NEW']
        self.sj_dict = {}
    
    def run(self):
        print(self.display_options())
        for sj_file in self.input:
            self.add_sj_file_to_dict(sj_file)
        
        self.write_sj_dict_to_file()
    
    def add_sj_file_to_dict(self, sj_file):
        filetype = None
        sjf = open(sj_file)
        for line in sjf:
            sj = SJobject(line, filetype)
            if filetype is None:
                filetype = sj.linetype
            
            if sj.hash in self.sj_dict:
                self.sj_dict[sj.hash].merge(sj)
            else:
                self.sj_dict[sj.hash] = sj
        
        sjf.close()
    
    def write_sj_dict_to_file(self):
        outfile_name = self.output
        if self.format == 'star' and not outfile_name.endswith('.SJ.out.tab'):
            outfile_name += '.SJ.out.tab'
        elif self.format == 'bed' and not outfile_name.endswith('.bed'):
            outfile_name += '.bed'

        output_file = open(outfile_name, 'w')
        for sj_hash, sj in sorted(list(self.sj_dict.items())):
            if self.new and not sj.new: continue
            if sj.count < self.min_reps: continue
            if sj.unique < self.min_unique: continue
            output_file.write(sj.as_string(self.format))
        
        output_file.close()
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend merge-sj |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        options_string += "  *** Parameters ***\n"
        options_string += "  Output format: {}\n".format(self.format)
        options_string += "  Min unique:    {}\n".format(self.min_unique)
        options_string += "  Min reps:      {}\n".format(self.min_reps)
        options_string += "  Only new SJs:  {}\n".format(self.new)
        return options_string

if __name__ == '__main__':
    from argument_parsers import sj_merge_parser as parser
    args = vars(parser.parse_args())
    obj = SJmerger(args)
    sys.exit(obj.run())
