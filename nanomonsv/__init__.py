#! /usr/bin/env python

from .arg_parser import create_parser

def main():

    parser = create_parser()
    args = parser.parse_args()
    args.func(args)

