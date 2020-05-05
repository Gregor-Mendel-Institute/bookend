#!/usr/bin/env python
# coding: utf-8
"""
    
    /|| bookend ||\\
    ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
    End-guided transcript assembly
    for short and long read RNA-seq
"""

import argparse
import logging, logging.config
import os
import sys

from bookend.core.argument_parsers import main_parser, Helper
from bookend.__init__ import __version__,__updated__,__date__

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'default': {
            'format': '%(asctime)s %(levelname)s %(name)s %(message)s'
        },
    },
    'handlers': {
        'stdout':{
            'class' : 'logging.StreamHandler',
            'stream'  : 'ext://sys.stdout',
            'formatter': 'default',
        },
        'stderr':{
            'class' : 'logging.StreamHandler',
            'stream'  : 'ext://sys.stderr',
            'level':'ERROR',
            'formatter': 'default',
        },
    },
    'root': {
        'handlers': ['stdout','stderr'],
        'level': 'INFO',
    },
}

logging.config.dictConfig(LOGGING)
log = logging.getLogger()

def get_parser(program_version_message):
    parser = main_parser
    parser.add_argument('-v', '--version', action='version', version=program_version_message)
    parser.add_argument('-l', '--loglevel',dest='log_level', help='Set log level', default='INFO', choices=['DEBUG','INFO','WARNING','ERROR'])
    return parser

def import_object(object_name):
    """Imports only the object needed to execute the subcommand."""
    if object_name == 'Assembler':
        from bookend.core.assemble import Assembler
        objectClass = Assembler
    elif object_name == 'Helper':
        from bookend.core.argument_parsers import Helper
        objectClass = Helper
    elif object_name == 'EndLabeler':
        from bookend.core.fastq_end_label import EndLabeler
        objectClass = EndLabeler
    elif object_name == 'EndLabeler':
        from bookend.core.bam_to_elr import BAMtoELRconverter
        objectClass = BAMtoELRconverter
    elif object_name == 'EndLabeler':
        from bookend.core.fasta_index import Indexer
        objectClass = Indexer
    elif object_name == 'EndLabeler':
        from bookend.core.bed_to_elr import BEDtoELRconverter
        objectClass = BEDtoELRconverter
    elif object_name == 'EndLabeler':
        from bookend.core.elr_to_bed import ELRtoBEDconverter
        objectClass = ELRtoBEDconverter
    elif object_name == 'EndLabeler':
        from bookend.core.gtf_merge import AnnotationMerger
        objectClass = AnnotationMerger
    elif object_name == 'EndLabeler':
        from bookend.core.elr_combine import ELRcombiner
        objectClass = ELRcombiner
    elif object_name == 'EndLabeler':
        from bookend.core.sam_sj_out import SAMtoSJconverter
        objectClass = SAMtoSJconverter
    elif object_name == 'EndLabeler':
        from bookend.core.sj_merge import SJmerger
        objectClass = SJmerger
    elif object_name == 'EndLabeler':
        from bookend.core.sj_to_bed import SJtoBEDconverter
        objectClass = SJtoBEDconverter
    elif object_name == 'EndLabeler':
        from bookend.core.elr_sort import ELRsorter
        objectClass = ELRsorter
    elif object_name == 'EndLabeler':
        from bookend.core.gtf_to_bed import GTFconverter
        objectClass = GTFconverter
    else:
        return None
    
    return objectClass

def main():
    """Passes commandline options to subfunctions."""
    version = 'v{}'.format(__version__)
    date = str(__updated__)
    program_version_message = '{} ({})'.format(version, date)
    
    parser = get_parser(program_version_message)    
    args = vars(parser.parse_args())
    try:
        log_level = args['log_level']
        log.setLevel(log_level)
        objectClass = import_object(args['object'])
        if objectClass is None:
            log.error('Subcommand not recognized. See bookend --help')
            return 1
        
        obj = objectClass(args)
        return obj.run()
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        log.exception(e)
        return 2

if __name__ == '__main__':
    sys.exit(main())
