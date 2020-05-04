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

from bookend.core.argument_parsers import main_parser
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
        obj = args['object'](args)
        return obj.run()
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        log.exception(e)
        return 2

if __name__ == '__main__':
    sys.exit(main())
