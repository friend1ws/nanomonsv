#! /usr/bin/env python

import pysam
import sys, math, re

def get_seq(reference, chr, start, end):

    seq = pysam.FastaFile(reference).fetch(chr, start - 1, end)
    """
    seq = ""    
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        if item[0] == ">": continue
        seq = seq + item.rstrip('\n').upper()
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')
    """

    if re.search(r'[^ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]', seq) is not None:
        print("The return value in get_seq function includes non-nucleotide characters:", file = sys.stderr)
        print(seq, file = sys.stderr)
        sys.exit(1)

    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


