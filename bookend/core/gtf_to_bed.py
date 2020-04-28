#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess

class GTFconverter:
    def __init__(self, args):
        """Converts each transcript of GTF-formatted input to BED12"""
        self.input = args['INPUT']
        self.output = args['OUT']
        
        self.command_string = 'bin/gtf-to-bed {}'.format(self.input)
        if self.output != 'stdout':
            self.command_string += ' > {}'.format(self.output)

    def run(self):
        if self.output != 'stdout':
            print(self.display_options())
        
        return subprocess.call(self.command_string, shell=True)
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend gtf-to-bed |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        return options_string

if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import gtf_to_bed_parser as parser
    args = vars(parser.parse_args())
    obj = GTFconverter(args)
    sys.exit(obj.run())
