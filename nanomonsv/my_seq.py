#! /usr/bin/env python

import pysam
import sys, math, re

from .logger import get_logger
logger = get_logger(__name__)

def get_seq(reference, chr, start, end):

    seq = pysam.FastaFile(reference).fetch(chr, start - 1, end)

    if re.search(r'[^ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]', seq) is not None:
        logger.error("The return value in get_seq function includes non-nucleotide characters: " % seq)
        sys.exit(1)

    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))
