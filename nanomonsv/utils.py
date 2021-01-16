#! /usr/bin/env python3

import sys, os
from subprocess import Popen, PIPE
import pysam

from .logger import get_logger
logger = get_logger(__name__)

def is_exists_bam(input_file):

    if input_file.startswith("s3://"):
        is_exists_s3(input_file)
    else:
        is_exists(input_file)


def is_exists(input_file):
    
    if not os.path.exists(input_file):
        logger.error("Input not exists: %s" % input_file)
        sys.exit(1)


def is_exists_s3(bam_object):

    from urllib.parse import urlparse
    import boto3

    obj_p = urlparse(bam_object)

    tbucket = obj_p.netloc
    tkey = obj_p.path
    if tkey.startswith("/"): tkey = tkey[1:]

    client = boto3.client("s3")
    try:
        response = client.head_object(Bucket = tbucket, Key = tkey)
    except:
        logger.error("Input not exists: %s: " % bam_object)
        sys.exit(1)


def is_exists_parsed_files(input_prefix):

    missing_files = []
    if not os.path.exists(input_prefix + ".deletion.sorted.bed.gz"):
        missing_files.append(input_prefix + ".deletion.sorted.bed.gz")
    if not os.path.exists(input_prefix + ".deletion.sorted.bed.gz.tbi"):
        missing_files.append(input_prefix + ".deletion.sorted.bed.gz.tbi")
    if not os.path.exists(input_prefix + ".insertion.sorted.bed.gz"):
        missing_files.append(input_prefix + ".insertion.sorted.bed.gz")
    if not os.path.exists(input_prefix + ".insertion.sorted.bed.gz.tbi"):
        missing_files.append(input_prefix + ".insertion.sorted.bed.gz.tbi")
    if not os.path.exists(input_prefix + ".rearrangement.sorted.bedpe.gz"):
        missing_files.append(input_prefix + ".rearrangement.sorted.bedpe.gz")
    if not os.path.exists(input_prefix + ".rearrangement.sorted.bedpe.gz.tbi"):
        missing_files.append(input_prefix + ".rearrangement.sorted.bedpe.gz.tbi")

    if len(missing_files) > 0:
        logger.error("One or more files via the nanomonsv parse stage do not exist. " + \
                     "Please check the prefix path. Also confirm whether nanomonsv parse ended properly. " +  \
                     "Missing files: " + ' '.join(missing_files))
        sys.exit(1)


def bam_format_check(bam_file):

    try:
        bam_t = pysam.AlignmentFile(bam_file)
    except:
        logger.error("BAM format error: %s" % bam_file)
        sys.exit(1)

def fasta_format_check(fasta_file):

    try:
        fasta_t = pysam.FastaFile(fasta_file)
    except:
        logger.error("FASTA format error: %s" % fasta_file)
        sys.exit(1)
       
 
def is_tool(executable):

    from shutil import which
    if which(executable) is None:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True


def libssw_check():

    # check whether libssw.so is in LD_LIBRARY_PATH
    sLibPath = ""
    for ld_path in os.environ["LD_LIBRARY_PATH"].split(':'):
        # print ld_path
        if os.path.exists(ld_path + "/libssw.so"):
            sLibPath = ld_path # + "/libssw.so"
            break
    if sLibPath == "":
        logger.error("Cannot find libssw.so in LD_LIBRARY_PATH")
        sys.exit(1)

