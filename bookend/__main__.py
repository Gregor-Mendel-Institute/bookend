#!/usr/bin/env python
# coding: utf-8
"""
    
    /|| bookend ||\\
    ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
    End-guided transcript assembly
    for short and long read RNA-seq
"""

import argparse
import os
import sys

from .core.argument_parsers import main_parser as parser
from .__init__ import __version__,__updated__,__date__

def get_parser(program_version_message):
    parser.add_argument('-v', '--version', action='version', version=program_version_message)
    return parser

def import_object(object_name):
    """Imports only the object needed to execute the subcommand."""
    if object_name == 'Assembler':
        from .core.assemble import Assembler
        objectClass = Assembler
    elif object_name == 'Helper':
        from .core.argument_parsers import Helper
        objectClass = Helper
    elif object_name == 'EndLabeler':
        from .core.fastq_end_label import EndLabeler
        objectClass = EndLabeler
    elif object_name == 'Condenser':
        from .core.elr_condense import Condenser
        objectClass = Condenser
    elif object_name == 'BAMtoELRconverter':
        from .core.bam_to_elr import BAMtoELRconverter
        objectClass = BAMtoELRconverter
    elif object_name == 'Indexer':
        from .core.fasta_index import Indexer
        objectClass = Indexer
    elif object_name == 'BEDtoELRconverter':
        from .core.bed_to_elr import BEDtoELRconverter
        objectClass = BEDtoELRconverter
    elif object_name == 'ELRtoBEDconverter':
        from .core.elr_to_bed import ELRtoBEDconverter
        objectClass = ELRtoBEDconverter
    elif object_name == 'AnnotationMerger':
        from .core.gtf_merge import AnnotationMerger
        objectClass = AnnotationMerger
    elif object_name == 'AssemblyClassifier':
        from .core.gtf_classify import AssemblyClassifier
        objectClass = AssemblyClassifier
    elif object_name == 'ELRcombiner':
        from .core.elr_combine import ELRcombiner
        objectClass = ELRcombiner
    elif object_name == 'SAMtoSJconverter':
        from .core.sam_sj_out import SAMtoSJconverter
        objectClass = SAMtoSJconverter
    elif object_name == 'SJmerger':
        from .core.sj_merge import SJmerger
        objectClass = SJmerger
    elif object_name == 'SJtoBEDconverter':
        from .core.sj_to_bed import SJtoBEDconverter
        objectClass = SJtoBEDconverter
    elif object_name == 'ELRsorter':
        from .core.elr_sort import ELRsorter
        objectClass = ELRsorter
    elif object_name == 'GTFconverter':
        from .core.gtf_to_bed import GTFconverter
        objectClass = GTFconverter
    elif object_name == 'ELRsimulator':
        from .core.elr_simulate import ELRsimulator
        objectClass = ELRsimulator
    else:
        from .core.argument_parsers import Helper
        objectClass = Helper
    
    return objectClass

def main():
    """Passes commandline options to subfunctions."""
    version = 'v{}'.format(__version__)
    date = str(__updated__)
    program_version_message = '{} ({})'.format(version, date)
    
    parser = get_parser(program_version_message)    
    args = vars(parser.parse_args())
    try:
        objectClass = import_object(args['object'])
        if objectClass is None:
            print('Subcommand not recognized. See bookend --help')
            return 1
        
        obj = objectClass(args)
        return obj.run()
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        print(e)
        return 2

if __name__ == '__main__':
    main()
