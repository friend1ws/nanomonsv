#! /usr/bin/env python

import subprocess
from .parse import *
from .filt import *
from .identify import *
from .long_read_validate import *

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

    cluster_rearrangement(args.tumor_prefix + ".rearrangement.sorted.bedpe.gz", args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe",
                     args.cluster_margin_size)

    filt_clustered_rearrangement1(args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe", args.tumor_prefix + ".rearrangement.sorted.clustered.filt1.bedpe",
                             args.min_tumor_variant_read_num, args.median_mapQ_thres, args.max_overhang_size_thres)

    filt_clustered_rearrangement2(args.tumor_prefix + ".rearrangement.sorted.clustered.filt1.bedpe", args.tumor_prefix + ".rearrangement.sorted.clustered.filt2.bedpe", 
                             args.control_prefix + ".rearrangement.sorted.bedpe.gz")


    cluster_insertion_deletion(args.tumor_prefix + ".deletion.sorted.bed.gz", args.tumor_prefix + ".deletion.sorted.clustered.bedpe")

    filt_clustered_insertion_deletion1(args.tumor_prefix + ".deletion.sorted.clustered.bedpe", args.tumor_prefix + ".deletion.sorted.clustered.filt1.bedpe",
                                       args.min_tumor_variant_read_num, args.median_mapQ_thres, args.max_overhang_size_thres)

    filt_clustered_insertion_deletion2(args.tumor_prefix + ".deletion.sorted.clustered.filt1.bedpe", args.tumor_prefix + ".deletion.sorted.clustered.filt2.bedpe",
                                       args.control_prefix + ".deletion.sorted.bed.gz")

    cluster_insertion_deletion(args.tumor_prefix + ".insertion.sorted.bed.gz", args.tumor_prefix + ".insertion.sorted.clustered.bedpe")

    filt_clustered_insertion_deletion1(args.tumor_prefix + ".insertion.sorted.clustered.bedpe", args.tumor_prefix + ".insertion.sorted.clustered.filt1.bedpe",
                                       args.min_tumor_variant_read_num, args.median_mapQ_thres, args.max_overhang_size_thres)

    filt_clustered_insertion_deletion2(args.tumor_prefix + ".insertion.sorted.clustered.filt1.bedpe", args.tumor_prefix + ".insertion.sorted.clustered.filt2.bedpe",
                                       args.control_prefix + ".insertion.sorted.bed.gz")

    identify(args.tumor_prefix + ".rearrangement.sorted.clustered.filt2.bedpe", 
             args.tumor_prefix + ".insertion.sorted.clustered.filt2.bedpe",
             args.tumor_prefix + ".deletion.sorted.clustered.filt2.bedpe",
             args.tumor_prefix + ".refined_bp.txt", args.tumor_bam, args.reference_fasta, args.debug)

    long_read_validate_main(args.tumor_prefix + ".refined_bp.txt",
                  args.tumor_bam,
                  args.tumor_prefix + ".refined_bp.validated.txt",
                  args.tumor_prefix + ".validated.tumor_sread.txt",
                  args.reference_fasta, 
                  args.control_bam, args.debug)

    is_control = True if args.control_bam is not None else False

    filt_final(args.tumor_prefix + ".refined_bp.validated.txt",
               args.tumor_prefix + ".validated.tumor_sread.txt",
               args.tumor_prefix + ".nanomonsv.result.txt",
               args.tumor_prefix + ".nanomonsv.supporting_read.txt",
               args.min_tumor_variant_read_num, args.min_tumor_VAF, args.max_control_variant_read_num, args.max_control_VAF, is_control)

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".rearrangement.sorted.clustered.filt1.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".rearrangement.sorted.clustered.filt2.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".deletion.sorted.clustered.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".deletion.sorted.clustered.filt1.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".deletion.sorted.clustered.filt2.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".insertion.sorted.clustered.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".insertion.sorted.clustered.filt1.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".insertion.sorted.clustered.filt2.bedpe"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".refined_bp.txt"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".refined_bp.validated.txt"])
        subprocess.check_call(["rm", "-rf", args.tumor_prefix + ".validated.tumor_sread.txt"])


def validate_main(args):

    long_read_validate_main(args.sv_list,
                            args.tumor_bam,
                            args.output + ".validated.txt",
                            args.output + ".validated.tumor_sread.txt",
                            args.reference_fasta,
                            args.control_bam, 
                            args.debug)

    is_control = True if args.control_bam is not None else False

    filt_final(args.output + ".validated.txt",
               args.output + ".validated.tumor_sread.txt",
               args.output, 
               args.output + ".supporting_read.txt",
               0, 0, float("inf"), float("inf"), is_control)

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output + ".validated.txt"])
        subprocess.check_call(["rm", "-rf", args.output + ".validated.tumor_sread.txt"])
    
