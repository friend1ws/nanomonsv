#! /usr/bin/env python

import subprocess
from .parse import *

def parse_main(args):

    parse_alignment_info(args.bam_file, args.output_prefix + ".tmp.alignment_info.txt")

    hout = open(args.output_prefix + ".tmp.alignment_info.name_sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", args.output_prefix + ".tmp.alignment_info.txt"], stdout = hout)
    hout.close()

    extract_bedpe_junction(args.output_prefix + ".tmp.alignment_info.name_sorted.txt", 
                           args.output_prefix + ".tmp.junction.bedpe",
                           args.split_alignment_check_margin, args.minimum_breakpoint_ambiguity)

    hout = open(args.output_prefix + ".tmp.junction.sorted.bedpe", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", 
                           args.output_prefix + ".tmp.junction.bedpe"], stdout = hout)
    hout.close()
  
    hout = open(args.output_prefix + ".junction.sorted.bedpe.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".tmp.junction.sorted.bedpe"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".junction.sorted.bedpe.gz"])


    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.alignment_info.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.alignment_info.name_sorted.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.junction.bedpe"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".junction.sorted.bedpe.gz"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.junction.sorted.bedpe"])

