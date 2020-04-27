#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess

class ELRsorter:
    def __init__(self, args):
        """Converts each line of BED-formatted input to ELR"""
        self.input = args['INPUT']
        self.output = args['OUT']
        self.ram = args['RAM']
        self.cpus = args['CPUS']
        self.temp = args['TEMP']

        self.command_string = 'bin/elr-sort --ram {} --cpus {} -T {} {}'.format(self.ram, self.cpus, self.temp, self.input)
        if self.output != 'stdout':
            self.command_string += ' > {}'.format(self.output)

    def run(self):
        if self.output != 'stdout':
            print(self.display_options())
        
        return subprocess.call(self.command_string, shell=True)
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend bed-to-elr |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Input file:    {}\n".format(self.input)
        options_string += "  Output file:   {}\n".format(self.output)
        options_string += "  *** Parameters ***\n"
        options_string += "  Temp directory: {}\n".format(self.temp)
        options_string += "  --ram: {}\n".format(self.ram)
        options_string += "  --cpus: {}\n".format(self.cpus)
        return options_string

if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import elr_sort_parser as parser
    args = vars(parser.parse_args())
    obj = ELRsorter(args)
    sys.exit(obj.run())
