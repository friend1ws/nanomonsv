#! /usr/bin/env python

import subprocess
from .parse import *
from .filt import *

def parse_main(args):

    parse_alignment_info(args.bam_file, args.output_prefix + ".tmp.deletion_info.txt", 
                                        args.output_prefix + ".tmp.insertion_info.txt", 
                                        args.output_prefix + ".tmp.rearrangement_info.txt")

    ####################
    # deletion processing
    hout = open(args.output_prefix + ".tmp.deletion.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", args.output_prefix + ".tmp.deletion_info.txt"], stdout = hout)
    hout.close()

    hout = open(args.output_prefix + ".deletion.sorted.bed.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".tmp.deletion.sorted.bed"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".deletion.sorted.bed.gz"])
    ####################

    ####################
    # insertion processing
    hout = open(args.output_prefix + ".tmp.insertion.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", args.output_prefix + ".tmp.insertion_info.txt"], stdout = hout)
    hout.close()
    
    hout = open(args.output_prefix + ".insertion.sorted.bed.gz", 'w') 
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".tmp.insertion.sorted.bed"], stdout = hout)
    hout.close()
     
    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".insertion.sorted.bed.gz"])
    ####################

    ####################
    # rearrangement processing
    hout = open(args.output_prefix + ".tmp.rearrangement_info.name_sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", args.output_prefix + ".tmp.rearrangement_info.txt"], stdout = hout)
    hout.close()

    extract_bedpe_junction(args.output_prefix + ".tmp.rearrangement_info.name_sorted.txt", 
                           args.output_prefix + ".tmp.rearrangement.bedpe",
                           args.split_alignment_check_margin, args.minimum_breakpoint_ambiguity)

    hout = open(args.output_prefix + ".tmp.rearrangement.sorted.bedpe", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", 
                           args.output_prefix + ".tmp.rearrangement.bedpe"], stdout = hout)
    hout.close()
  
    hout = open(args.output_prefix + ".rearrangement.sorted.bedpe.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".tmp.rearrangement.sorted.bedpe"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".rearrangement.sorted.bedpe.gz"])
    ####################

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.deletion_info.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.deletion.sorted.bed"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.insertion_info.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.insertion.sorted.bed"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.rearrangement_info.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.rearrangement_info.name_sorted.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.rearrangement.bedpe"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.rearrangement.sorted.bedpe"])


def get_main(args):

    cluster_junction(args.tumor_prefix + ".junction.sorted.bedpe.gz", args.tumor_prefix + ".junction.sorted.clustered.bedpe",
                     args.cluster_margin_size)

    filt_clustered_junction1(args.tumor_prefix + ".junction.sorted.clustered.bedpe", args.tumor_prefix + ".junction.sorted.clustered.filt1.bedpe",
                             args.read_num_thres, args.median_mapQ_thres, args.max_overhang_size_thres)

    filt_clustered_junction2(args.tumor_prefix + ".junction.sorted.clustered.filt1.bedpe", args.tumor_prefix + ".junction.sorted.clustered.filt2.bedpe", 
                             args.control_prefix + ".junction.sorted.bedpe.gz")


