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
        logger.error("BAM format error!")
        sys.exit(1)

        
def is_tool(executable):

    from shutil import which
    if which(executable) is None:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True

