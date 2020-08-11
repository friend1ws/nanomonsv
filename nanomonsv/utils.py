#! /usr/bin/env python

import sys
from subprocess import Popen, PIPE
import pysam

from .logger import get_logger
logger = get_logger(__name__)

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

