#! /usr/bin/env python

import sys
from .arg_parser import create_parser

def main():

    parser = create_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args.func(args)

