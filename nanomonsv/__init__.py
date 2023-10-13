#! /usr/bin/env python

import sys
from .arg_parser import create_parser
from .logger import set_logger_dir

def main():

    parser = create_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    set_logger_dir(args.out_folder)
    args.func(args)
