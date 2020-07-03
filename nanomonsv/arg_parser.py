#! /usr/bin/env python

import argparse

from .run import *
from .version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "nanomonsv")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    subparsers = parser.add_subparsers()

    ##########
    # parse
    parse = subparsers.add_parser("parse",
                                help = "Parse supporting reads for candidate structural variations")


    parse.add_argument("bam_file", default = None, type = str,
                     help = "Path to input BAM file")

    parse.add_argument("output_prefix", type = str,
                       help = "Prefix of output files")

    parse.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    parse.add_argument("--split_alignment_check_margin", default = 50, type = int,
                       help = "Two split alignments whose margin sizes are no more than this value is counted as candidate breakpoint")

    parse.add_argument("--minimum_breakpoint_ambiguity", default = 20, type = int,
                       help = "Sizes of ambiguities of breakpoint positions from the observed ones") 

    parse.set_defaults(func = parse_main)
    ##########
    
    ##########
    # get
    get = subparsers.add_parser("get",
                                help = "List up reliable structural variations with refined breakpoint positions")

    get.add_argument("tumor_prefix", type = str,
                      help = "Prefix of tumor data processed in parse step")
       
    get.add_argument("tumor_bam", default = None, type = str,
                      help = "Path to tumor BAM file")
 
    get.add_argument("reference_fasta", metavar = "reference.fa", default = None, type = str,
                     help = "the path to the reference genome sequence")

    get.add_argument("--control_prefix", type = str,
                     help = "Prefix of matched control data processed in parse step")

    get.add_argument("--control_bam", type = str,
                     help = "Path to control BAM file")

    get.add_argument("--min_tumor_variant_read_num", default = 3, type = int,
                     help = "Minimum required supporting read number for a tumor sample")

    get.add_argument("--min_tumor_VAF", default = 0.05, type = float,
                     help = "Minimum required variant allele frequency for a tumor sample")

    get.add_argument("--max_control_variant_read_num", default = 1, type = int,
                     help = "Maximum allowed supporting read number for a control sample")
        
    get.add_argument("--max_control_VAF", default = 0.03, type = float,
                     help = "Maximum allowed variant allele frequeycy for a control sample")

    get.add_argument("--cluster_margin_size", default = 100, type = int,
                     help = "Two breakpoints are margined if they are within this threshould value")

    # get.add_argument("--read_num_thres", default = 3, type = int,
    #                  help = "Minimum required supporting read number for structural variation candidates")

    get.add_argument("--median_mapQ_thres", default = 20, type = int,
                     help = "Threshould for median mapping quality")

    get.add_argument("--max_overhang_size_thres", default = 100, type = int,
                     help = "Threshould for maximum overhang size")

    get.add_argument("--var_read_min_mapq", default = 0, type = int,
                     help = "Threshould for mapping quality in validate step")

    # get.add_argument("--control_read_num_thres", default = 0, type = int,
    #                  help = "Filter if the number of supporting reads for the control sample is larger than this value")

    get.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    get.set_defaults(func = get_main)
    ##########

    """
    ##########
    # validate 
    validate = subparsers.add_parser("validate",
                                     help = "Validate GenomonSV format SV list using nanopore sequence reads")
 
    validate.add_argument("sv_list_file", default = None, type = str,
                          help = "Path to GenomonSV format SV list file")

    validate.add_argument("tumor_bam", default = None, type = str,
                          help = "Path to tumor BAM file")

    validate.add_argument("output", type = str,
                       help = "Path to output file")

    validate.add_argument("reference_fasta", metavar = "reference.fa", default = None, type = str,
                          help = "Path to the reference genome sequence")

    validate.add_argument("--control_bam", default = None, type = str,
                          help = "Path to control BAM file")

    validate.add_argument("--var_read_min_mapq", default = 40, type = int,
                          help = "Threshould for mapping quality in validate step")

    validate.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    validate.set_defaults(func = validate_main)
    ##########
    """

    ##########
    # insert_classify
    insert_classify = subparsers.add_parser("insert_classify",
                                            help = "Classify long insertion into LINE1, Alu, SVA, and so on")
    
    insert_classify.add_argument("sv_list_file", default = None, type =str,
                                 help = "Path to nanomonsv get result file")
    
    insert_classify.add_argument("output_file", type = str,
                                 help = "Path to output file")

    insert_classify.add_argument("reference_fasta", metavar = "reference.fa", default = None, type = str,
                                 help = "Path to the reference genome sequence")

    insert_classify.add_argument("--grc", default = False, action = 'store_true',
                                 help = "Deprecated. This is not used any more. Convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")

    insert_classify.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                                 help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    insert_classify.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    insert_classify.set_defaults(func = insert_classify_main)
    ##########

    return parser
  


    
