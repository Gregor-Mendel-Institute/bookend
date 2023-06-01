#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
if __name__ == '__main__':
    sys.path.append('../../bookend')

class SJtoBEDconverter:
    def __init__(self, args):
        """Converts an SJ.out.tab to an SJ.bed file."""
        self.input = args['INPUT']
        self.output = args['OUT']
    
    def run(self):
        counter=0
        output_file = open(self.output, 'w')
        for line in open(self.input):
            counter+=1
            l=line.rstrip().split('\t')
            reads = str(int(l[6])+int(l[7]))
            if l[3]=='1':
                outline = '\t'.join([l[0],str(int(l[1])-1),l[2],'SJ.'+str(counter),reads,'+'])
            elif l[3]=='2':
                outline = '\t'.join([l[0],str(int(l[1])-1),l[2],'SJ.'+str(counter),reads,'-'])

            output_file.write(outline+'\n')
        
        output_file.close()

if __name__ == '__main__':
    from argument_parsers import sj_to_bed_parser as parser
    args = vars(parser.parse_args())
    obj = SJtoBEDconverter(args)
    sys.exit(obj.run())
