#! /usr/bin/env python

import pysam
import math

def get_seq(reference, chr, start, end):

    seq = ""    
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        if item[0] == ">": continue
        seq = seq + item.rstrip('\n').upper()
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


