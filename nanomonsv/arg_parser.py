#! /usr/bin/env python3

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


    parse.add_argument("alignment_file", default = None, type = str,
                     help = "Path to input BAM or CRAM file")

    parse.add_argument("output_prefix", type = str,
                       help = "Prefix of output files")

    parse.add_argument("--reference_fasta", metavar = "reference.fa", default = None, type = str,
                     help = "the path to the reference genome sequence")

    parse.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    parse.add_argument("--split_alignment_check_margin", default = 50, type = int,
                       help = "Two split alignments whose margin sizes are no more than this value is counted as candidate breakpoint (default: 50)")

    parse.add_argument("--minimum_breakpoint_ambiguity", default = 20, type = int,
                       help = "Sizes of ambiguities of breakpoint positions from the observed ones (default: 20)") 

    parse.set_defaults(func = parse_main)
    ##########
   
    ##########
    # merge_control
    merge_control = subparsers.add_parser("merge_control", 
                                          help = "Merge control panel files from parse results")

    merge_control.add_argument("prefix_list_file", default = None, type = str,
                               help = "Path to prefix for the parse function for control panels")
        
    merge_control.add_argument("output_prefix", default = None, type = str,
                               help = "Prefix of merged control files")

    merge_control.set_defaults(func = merge_control_main)
    ##########

    ##########
    # get
    get = subparsers.add_parser("get",
                                help = "List up reliable structural variations with refined breakpoint positions")

    get.add_argument("tumor_prefix", type = str,
                      help = "Prefix of tumor data processed in parse step")
       
    get.add_argument("tumor_bam", default = None, type = str,
                      help = "Path to tumor alignment (BAM or CRAM) file")
 
    get.add_argument("reference_fasta", metavar = "reference.fa", default = None, type = str,
                     help = "the path to the reference genome sequence")

    get.add_argument("--control_prefix", type = str, # required = True,
                     help = "Prefix of matched control data processed in parse step")

    get.add_argument("--control_bam", type = str, # required = True,
                     help = "Path to control alignment (BAM or CRAM) file")

    get.add_argument("--control_panel_prefix", type = str, 
                     help = "Prefix of non-matched control panel data processed in merge_control step")

    get.add_argument("--min_tumor_variant_read_num", default = 3, type = int,
                     help = "Minimum required supporting read number for a tumor sample (default: 3)")

    get.add_argument("--min_tumor_VAF", default = 0.05, type = float,
                     help = "Minimum required variant allele frequency for a tumor sample (default: 0.05)")

    get.add_argument("--max_control_variant_read_num", default = 1, type = int,
                     help = "Maximum allowed supporting read number for a control sample (default: 1)")
        
    get.add_argument("--max_control_VAF", default = 0.03, type = float,
                     help = "Maximum allowed variant allele frequeycy for a control sample (default: 0.03)")

    get.add_argument("--min_indel_size", default = 50, type = int,
                     help = "Minimum indel size for the output (default: 50)")

    get.add_argument("--max_panel_read_num", default = 1, type = int,
                     help= "Maximum allowed supporting read number for a nonmatched control sample (default: 1)")

    get.add_argument("--max_panel_sample_num", default = 0, type = int,
                     help= "Maximum allowed sample number for a nonmatched control sample (default: 0)")

    get.add_argument("--cluster_margin_size", default = 100, type = int,
                     help = "Two breakpoints are margined if they are within this threshould value (default: 100)")

    get.add_argument("--median_mapQ_thres", default = 20, type = int,
                     help = "Threshould for median mapping quality (default: 20)")

    get.add_argument("--max_overhang_size_thres", default = 100, type = int,
                     help = "Threshould for maximum overhang size (default: 100)")

    get.add_argument("--check_read_max_num", default = 500, type = int,
                     help = "The maximum number of reads to check per breakpoint in the phase of realignment validation (default: 500)")

    get.add_argument("--var_read_min_mapq", default = 0, type = int,
                     help = "Threshould for mapping quality in validate step (default: 0)")

    get.add_argument("--validation_score_ratio_thres", default = 1.2, type = float,
                     help = "Threshould for threshould for SV segment validation by alignment")

    get.add_argument("--sw_jump_params", nargs = 4, default = [1, 3, 3, 2], type = int,
                     help = "Parameters (match score, mismatch penalty, gap penalty, insertion penalty) for one-time smith-waterman algorithm (default: [1, 3, 3, 2]")

    # get.add_argument("--use_ssw_lib", default = False, action = 'store_true',
    #                  help = "Use SSW Library. This is for backward comaptibility, and may be removed in the future (default: False)")

    get.add_argument("--qv10", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 10. Recommended for ONT data called by Guppy before version 5")

    get.add_argument("--qv15", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 15. Recommended for ONT data called by Guppy version 5, 6.")

    get.add_argument("--qv20", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 20. Recommended for ONT data with Q20+ chemistry.")

    get.add_argument("--qv25", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality above 25. Recommended for PacBio Hifi data.")

    get.add_argument("--use_racon", default = False, action = 'store_true',
                     help = "Use racon for error correction of clustered putative supporting reads (default: False)")

    get.add_argument("--single_bnd", default = False, action = 'store_true',
                     help = "Generate single end breakpoints (default: False)")

    # get.add_argument("--control_read_num_thres", default = 0, type = int,
    #                  help = "Filter if the number of supporting reads for the control sample is larger than this value")

    # get.add_argument("--threads", default = 1, type = int,
    #                  help = "Number of parallel threads to use (not recommended) (default: 1)")

    get.add_argument("--processes", default = 1, type = int,
                     help = "Number of parallel processes to use (default: 1)")

    get.add_argument("--sort_option", type = str, default = "-S 1G", 
                     help = "Options for Linux sort command (default: '-S 1G')")

    get.add_argument("--max_memory_minimap2", type = int, default = 2, 
                     help = "Maximum memory size (Gbyte) for minimap2 (default: 2)")

    get.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files (default: False)")

    get.set_defaults(func = get_main)
    ##########

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

    validate.add_argument("--validation_score_ratio_thres", default = 1.2, type = float, 
                     help = "Threshould for threshould for SV segment validation by alignment")

    validate.add_argument("--qv10", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 10. Recommended for ONT data called by Guppy before version 5")
    
    validate.add_argument("--qv15", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 15. Recommended for ONT data called by Guppy version 5, 6.")
    
    validate.add_argument("--qv20", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality of around 20. Recommended for ONT data with Q20+ chemistry.")

    validate.add_argument("--qv25", default = False, action = 'store_true',
                     help = "Parameter preset for sequencing data with a base quality above 25. Recommended for PacBio Hifi data.")

    validate.add_argument("--sort_option", metavar = "-S 1G", type = str, default = "-S 2G", 
                     help = "options for sort command")

    validate.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    validate.set_defaults(func = validate_main)
    ##########

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

    # insert_classify.add_argument("--grc", default = False, action = 'store_true',
    #                              help = "Deprecated. This is not used any more. Convert chromosome names to Genome Reference Consortium nomenclature (default: %(default)s)")
    insert_classify.add_argument("gtf_file", metavar = "gencode.gtf.gz", default = None, type = str,
                                 help = "Path to GFT file for transcript")

    # insert_classify.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg38",
    #                              help = "Genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    insert_classify.add_argument("LINE1_db", default = None, type = str,
                                 help = "Path to LINE1 position file")

    insert_classify.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    insert_classify.set_defaults(func = insert_classify_main)
    ##########
    
    ##########
    # connect
    connect = subparsers.add_parser("connect",
                                    help = "Connect called SVs by support reads")

    connect.add_argument("nanomonsv_result_file", default = None, type = str,
                         help = "Path to nanomonsv get result file")

    connect.add_argument("support_read_file", type = str,
                         help = "PATH to support read file")
    
    connect.add_argument("output_prefix", type = str,
                         help = "Prefix of output files")

    connect.add_argument("--consistent_pos_margin", default = 1.3, type = float,
                         help = "Consistency coefficient of distances of read and reference between two SVs")
    
    connect.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    connect.set_defaults(func = connect_main)
    ##########

    return parser


    
