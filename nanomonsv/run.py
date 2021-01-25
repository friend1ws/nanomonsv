#! /usr/bin/env python3

import os, subprocess, shutil
from .parse import *
from .cluster import *
from .cluster_sbnd import *
from .gather_support_read_seq import *
from .generate_consensus import *
from .locate_bp import *
from .locate_bp_sbnd import *
from .count_sread_by_alignment import *
from .post_proc import *
from .insert_classify import *
from .utils import *
from .logger import get_logger


logger = get_logger(__name__)

def parse_main(args):

    # check if the executables exist
    is_tool("tabix")
    is_tool("bgzip")

    # check input file existences
    is_exists_bam(args.bam_file)

    # BAM format check
    bam_format_check(args.bam_file)

    # make directory for the output prefix
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir != '' and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ####################
    parse_alignment_info(args.bam_file, args.output_prefix + ".tmp.deletion_info.txt", 
                                        args.output_prefix + ".tmp.insertion_info.txt", 
                                        args.output_prefix + ".tmp.rearrangement_info.txt",
                                        args.output_prefix + ".tmp.bp_info.txt")

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
                           args.output_prefix + ".tmp.rearrangement.bedpe") # ,
                           # args.split_alignment_check_margin, args.minimum_breakpoint_ambiguity)

    hout = open(args.output_prefix + ".tmp.rearrangement.sorted.bedpe", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", 
                           args.output_prefix + ".tmp.rearrangement.bedpe"], stdout = hout)
    hout.close()
  
    hout = open(args.output_prefix + ".rearrangement.sorted.bedpe.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".tmp.rearrangement.sorted.bedpe"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".rearrangement.sorted.bedpe.gz"])
    ####################

    ####################
    # breakpoint processing
    hout = open(args.output_prefix + ".bp_info.sorted.bed", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", args.output_prefix + ".tmp.bp_info.txt"], stdout = hout)
    hout.close()

    hout = open(args.output_prefix + ".bp_info.sorted.bed.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", args.output_prefix + ".bp_info.sorted.bed"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", args.output_prefix + ".bp_info.sorted.bed.gz"])
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
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.bp_info.txt"])
        subprocess.check_call(["rm", "-rf", args.output_prefix + ".bp_info.sorted.bed"])


def get_main(args):

    # check if the executables exist
    if args.use_racon: 
        is_tool("racon")
    else:
        is_tool("mafft")
    if args.use_ssw_lib: libssw_check()

    # check existences
    is_exists_bam(args.tumor_bam)
    is_exists(args.reference_fasta)
    if args.control_bam is not None: is_exists_bam(args.control_bam)
   
    # check parsed files existences
    is_exists_parsed_files(args.tumor_prefix)
    if args.control_prefix is not None: is_exists_parsed_files(args.control_prefix)
 
    # BAM format check
    bam_format_check(args.tumor_bam)
    fasta_format_check(args.reference_fasta)
    if args.control_bam is not None: bam_format_check(args.control_bam)

    control_rearrangement_bedpe = None
    control_deletion_bed = None
    control_insertion_bed = None
    if args.control_prefix is not None: 
        control_rearrangement_bedpe = args.control_prefix + ".rearrangement.sorted.bedpe.gz"
        control_deletion_bed = args.control_prefix + ".deletion.sorted.bed.gz"
        control_insertion_bed = args.control_prefix + ".insertion.sorted.bed.gz"
        if args.single_bnd is not None:
            control_bp_bed = args.control_prefix + ".bp_info.sorted.bed.gz"

    bp_bed = None
    if args.use_racon:
        bp_bed = args.tumor_prefix + ".bp_info.sorted.bed.gz"
    ####################
    logger.info("Clustering rearrangement type supporting reads for putative SVs")
    cluster_supporting_reads(args.tumor_prefix + ".rearrangement.sorted.bedpe.gz", 
        args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe",
        "rearrangement", control_junction_bedpe = control_rearrangement_bedpe,
        read_num_thres = args.min_tumor_variant_read_num, cluster_margin_size = args.cluster_margin_size, 
        median_mapQ_thres = args.median_mapQ_thres, max_overhang_size_thres = args.max_overhang_size_thres,
        maximum_local_variant_num = 1000, debug = args.debug)

    logger.info("Clustering insertion type supporting reads for putative SVs")
    cluster_supporting_reads(args.tumor_prefix + ".insertion.sorted.bed.gz",
        args.tumor_prefix + ".insertion.sorted.clustered.bedpe",
        "insertion", control_junction_bedpe = control_insertion_bed, bp_bed = bp_bed,
        read_num_thres = args.min_tumor_variant_read_num, cluster_margin_size = args.cluster_margin_size, 
        median_mapQ_thres = args.median_mapQ_thres, max_overhang_size_thres = args.max_overhang_size_thres,
        debug = args.debug)

    logger.info("Clustering deletion type supporting reads for putative SVs")
    cluster_supporting_reads(args.tumor_prefix + ".deletion.sorted.bed.gz",
        args.tumor_prefix + ".deletion.sorted.clustered.bedpe",
        "deletion", control_junction_bedpe = control_deletion_bed, bp_bed = bp_bed,
        read_num_thres = args.min_tumor_variant_read_num, cluster_margin_size = args.cluster_margin_size, 
        median_mapQ_thres = args.median_mapQ_thres, max_overhang_size_thres = args.max_overhang_size_thres,
        debug = args.debug)

    if args.single_bnd:
        logger.info("Clustering single breakend type supporting reads fro putative SVs")
        cluster_supporting_reads_sbnd(args.tumor_prefix + ".bp_info.sorted.bed.gz",
            args.tumor_prefix + ".singlebreakend.sorted.clustered.bed", control_bed = control_bp_bed,
            read_num_thres = args.min_tumor_variant_read_num, cluster_margin_size = args.cluster_margin_size,
            median_mapQ_thres = args.median_mapQ_thres, debug = args.debug)

    logger.info("Gathering sequences of supporting reads")
    if args.single_bnd:
        gather_support_read_seq(args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe",
            args.tumor_prefix + ".insertion.sorted.clustered.bedpe",
            args.tumor_prefix + ".deletion.sorted.clustered.bedpe",
            args.tumor_prefix + ".support_read_seq.txt",
            args.tumor_bam, single_breakend_file = args.tumor_prefix + ".singlebreakend.sorted.clustered.bed",
            output_file_sbind = args.tumor_prefix + ".support_read_seq.sbnd.txt" )
    else:
        gather_support_read_seq(args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe",
            args.tumor_prefix + ".insertion.sorted.clustered.bedpe",
            args.tumor_prefix + ".deletion.sorted.clustered.bedpe",
            args.tumor_prefix + ".support_read_seq.txt",
            args.tumor_bam)

    logger.info("Generate consensus sequences")
    generate_consensus(args.tumor_prefix + ".support_read_seq.txt",
        args.tumor_prefix + ".consensus_seq.txt",
        use_racon = args.use_racon, debug = args.debug)
    if args.single_bnd:
        generate_consensus_sbnd(args.tumor_prefix + ".support_read_seq.sbnd.txt",
            args.tumor_prefix + ".consensus_seq.sbnd.txt",
            use_racon = args.use_racon, debug = args.debug)
    
    logger.info("Locating single-base resolution break points for candidate SVs")
    locate_bp(args.tumor_prefix + ".consensus_seq.txt",
        args.tumor_prefix + ".refined_bp.txt",
        args.reference_fasta, args.debug)

    if args.single_bnd:
        locate_bp_sbnd(args.tumor_prefix + ".consensus_seq.sbnd.txt",
            args.tumor_prefix + ".refined_bp.sbnd.txt",
            args.reference_fasta, args.debug)
   
    logger.info("Counting the number of supprting read for the tumor by realignment of SV candidate segments")
    if args.single_bnd:
        count_sread_by_alignment(args.tumor_prefix + ".refined_bp.txt", args.tumor_bam, 
            args.tumor_prefix + ".realignment.tumor.sread_count.txt", args.tumor_prefix + ".realignment.tumor.sread_info.txt", args.reference_fasta,
            sbnd_file = args.tumor_prefix + ".refined_bp.sbnd.txt", output_count_file_sbnd = args.tumor_prefix + ".realignment.tumor.sread_count.sbnd.txt",
            output_alignment_info_file_sbnd = args.tumor_prefix + ".realignment.tumor.sread_info.sbnd.txt",
            var_read_min_mapq = args.var_read_min_mapq, use_ssw_lib = args.use_ssw_lib, debug = args.debug)
    else:
        count_sread_by_alignment(args.tumor_prefix + ".refined_bp.txt", args.tumor_bam,
            args.tumor_prefix + ".realignment.tumor.sread_count.txt", args.tumor_prefix + ".realignment.tumor.sread_info.txt",
            args.reference_fasta, args.var_read_min_mapq, args.use_ssw_lib, args.debug)
 
    if args.control_bam is not None:
        logger.info("Counting the number of supprting read for the control by realignment of SV candidate segments")
        if args.single_bnd:
            count_sread_by_alignment(args.tumor_prefix + ".refined_bp.txt", args.control_bam, 
                args.tumor_prefix + ".realignment.control.sread_count.txt", args.tumor_prefix + ".realignment.control.sread_info.txt", args.reference_fasta,
                sbnd_file = args.tumor_prefix + ".refined_bp.sbnd.txt", output_count_file_sbnd = args.tumor_prefix + ".realignment.control.sread_count.sbnd.txt",
                output_alignment_info_file_sbnd = args.tumor_prefix + ".realignment.control.sread_info.sbnd.txt",
                var_read_min_mapq = args.var_read_min_mapq, use_ssw_lib = args.use_ssw_lib, debug = args.debug)

        else:
            count_sread_by_alignment(args.tumor_prefix + ".refined_bp.txt", args.control_bam,
                args.tumor_prefix + ".realignment.control.sread_count.txt", args.tumor_prefix + ".realignment.control.sread_info.txt",
                args.reference_fasta, args.var_read_min_mapq, args.use_ssw_lib, args.debug)

    logger.info("Final processing") 
    control_sread_count_file = args.tumor_prefix + ".realignment.control.sread_count.txt" if args.control_bam is not None else None
    integrate_realignment_result(args.tumor_prefix + ".realignment.tumor.sread_count.txt", control_sread_count_file,
        args.tumor_prefix + ".nanomonsv.result.txt",
        args.min_tumor_variant_read_num, args.min_tumor_VAF, args.max_control_variant_read_num, args.max_control_VAF)

    proc_sread_info_file(args.tumor_prefix + ".realignment.tumor.sread_info.txt",
        args.tumor_prefix + ".nanomonsv.result.txt",
        args.tumor_prefix + ".nanomonsv.supporting_read.txt")
 
    if not args.debug:
        os.remove(args.tumor_prefix + ".rearrangement.sorted.clustered.bedpe")
        os.remove(args.tumor_prefix + ".insertion.sorted.clustered.bedpe")
        os.remove(args.tumor_prefix + ".deletion.sorted.clustered.bedpe")
        os.remove(args.tumor_prefix + ".support_read_seq.txt")
        os.remove(args.tumor_prefix + ".consensus_seq.txt")
        os.remove(args.tumor_prefix + ".refined_bp.txt")
        os.remove(args.tumor_prefix + ".realignment.tumor.sread_count.txt")
        os.remove(args.tumor_prefix + ".realignment.tumor.sread_info.txt")
        
        if args.control_bam is not None:
            os.remove(args.tumor_prefix + ".realignment.control.sread_count.txt")
            os.remove(args.tumor_prefix + ".realignment.control.sread_info.txt")


def validate_main(args):
   
    # executable check
    if args.use_ssw_lib: libssw_check()

    
    logger.info("Counting the number of supprting read for the tumor by realignment of SV candidate segments")
    count_sread_by_alignment(args.sv_list_file, args.tumor_bam,
        args.output + ".realignment.tumor.sread_count.txt", args.output + ".realignment.tumor.sread_info.txt",
        args.reference_fasta, args.var_read_min_mapq, args.use_ssw_lib, args.debug)

    if args.control_bam is not None:
        logger.info("Counting the number of supprting read for the control by realignment of SV candidate segments")
        count_sread_by_alignment(args.sv_list_file, args.control_bam,
            args.output + ".realignment.control.sread_count.txt", args.output + ".realignment.control.sread_info.txt",
            args.reference_fasta, args.var_read_min_mapq, args.use_ssw_lib, args.debug)

    logger.info("Final processing")
    control_sread_count_file = args.output + ".realignment.control.sread_count.txt" if args.control_bam is not None else None
    integrate_realignment_result(args.output + ".realignment.tumor.sread_count.txt", control_sread_count_file, args.output,
        0, 0, float("inf"), float("inf"))

    if not args.debug:
        os.remove(args.output + ".realignment.tumor.sread_count.txt")
        os.remove(args.output + ".realignment.tumor.sread_info.txt")
        os.remove(args.output + ".realignment.control.sread_count.txt")
        os.remove(args.output + ".realignment.control.sread_info.txt")
    ####################
    """
    long_read_validate_main(args.sv_list_file,
                            args.tumor_bam,
                            args.output + ".validated.txt",
                            args.output + ".validated.tumor_sread.txt",
                            args.reference_fasta,
                            args.control_bam, 
                            args.var_read_min_mapq,
                            args.use_ssw_lib, args.debug)

    is_control = True if args.control_bam is not None else False

    filt_final(args.output + ".validated.txt",
               args.output + ".validated.tumor_sread.txt",
               args.output, 
               args.output + ".supporting_read.txt",
               0, 0, float("inf"), float("inf"), True, is_control)

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output + ".validated.txt"])
        subprocess.check_call(["rm", "-rf", args.output + ".validated.tumor_sread.txt"])
    """
    

def insert_classify_main(args):

    import tempfile
    import annot_utils.exon

    # check if the executables exist
    is_tool("minimap2")
    is_tool("bedtools")
    is_tool("bwa")
    is_tool("RepeatMasker")

    make_fasta_file(args.sv_list_file, args.output_file + ".tmp.fasta", args.output_file + ".tmp.seq_id.txt")
   
    ##########
    # processed pseudo gene
    annot_utils.exon.make_exon_info(args.output_file + ".tmp.exon.bed.gz", "gencode", args.genome_id, args.grc, True)

    with open(args.output_file + ".tmp.minimap2.sam", 'w') as hout:
        subprocess.check_call(["minimap2", "-ax", "splice", args.reference_fasta, args.output_file + ".tmp.fasta"], stdout = hout)

    sam2bed_split(args.output_file + ".tmp.minimap2.sam", args.output_file + ".tmp.minimap2.filt.bed")
    
    with open(args.output_file + ".tmp.minimap2.filt.exon.bed", 'w') as hout:
        subprocess.check_call(["bedtools", "intersect", "-a", args.output_file + ".tmp.minimap2.filt.bed", 
                               "-b", args.output_file + ".tmp.exon.bed.gz", "-wo"], stdout = hout)
 
    pp_proc_filt_exon(args.output_file + ".tmp.minimap2.filt.exon.bed", 
                      args.output_file + ".tmp.seq_id.txt", 
                      args.output_file + ".tmp.ppseudo.txt")
    ##########

    ##########
    # repeat masker
    output_dir = os.path.dirname(args.output_file)
    tmpdir_rmsk = tempfile.mkdtemp()
    # tmpdir_rmsk = "tmp"
    subprocess.check_call(["RepeatMasker", "-species", "human", args.output_file + ".tmp.fasta", "-dir", tmpdir_rmsk])

    summarize_rmsk(tmpdir_rmsk + '/' + os.path.basename(args.output_file + ".tmp.fasta") + ".out", args.output_file + ".tmp.rmsk.txt")

    check_tsd_polyAT(args.output_file + ".tmp.fasta", args.output_file + ".tmp.seq_id.txt", 
                     args.reference_fasta, args.output_file + ".tmp.tsd.polyAT.txt")

    shutil.rmtree(tmpdir_rmsk)
    ##########

    ##########
    # alignment to reference genome
    with open(args.output_file + ".tmp.bwa.sam", 'w') as hout:
        print(' '.join(["bwa", "mem", "-h", "200", args.reference_fasta, args.output_file + ".tmp.fasta"]))
        subprocess.check_call(["bwa", "mem", "-h", "200", args.reference_fasta, args.output_file + ".tmp.fasta"], stdout = hout)

    summarize_bwa_alignment2(args.output_file + ".tmp.bwa.sam", args.output_file + ".tmp.seq_id.txt", args.output_file + ".tmp.alignment.txt")

    organize_info(args.output_file + ".tmp.rmsk.txt", args.output_file + ".tmp.alignment.txt", 
                 args.output_file + ".tmp.tsd.polyAT.txt", args.output_file + ".tmp.seq_id.txt", 
                 args.output_file + ".tmp.org.txt", args.genome_id)

    annotate_sv_file(args.sv_list_file, args.output_file + ".tmp.org.txt", args.output_file + ".tmp.ppseudo.txt",
                     args.output_file + ".tmp.seq_id.txt", args.output_file)

    if not args.debug:
        os.remove(args.output_file + ".tmp.fasta")
        os.remove(args.output_file + ".tmp.seq_id.txt")
        os.remove(args.output_file + ".tmp.bwa.sam")
        os.remove(args.output_file + ".tmp.exon.bed.gz")
        os.remove(args.output_file + ".tmp.exon.bed.gz.tbi")
        os.remove(args.output_file + ".tmp.minimap2.filt.bed")
        os.remove(args.output_file + ".tmp.minimap2.filt.exon.bed")
        os.remove(args.output_file + ".tmp.ppseudo.txt")
        os.remove(args.output_file + ".tmp.rmsk.txt")
        os.remove(args.output_file + ".tmp.tsd.polyAT.txt")
        os.remove(args.output_file + ".tmp.minimap2.sam")
        os.remove(args.output_file + ".tmp.alignment.txt")
        os.remove(args.output_file + ".tmp.org.txt")


