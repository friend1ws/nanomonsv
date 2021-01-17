#! /usr/bin/env python3

import sys, os, subprocess, statistics 
import pysam

from .logger import get_logger
from .my_seq import reverse_complement

logger = get_logger(__name__)

 
def set_readid2alignment(readid2alignment, input_file, mode, alignment_margin):

    # readid2alignment = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join([F[0], F[1], F[2], F[8], F[3], F[4], F[5], F[9], mode])
            readids = F[6].split(';')

            if mode == "r":
                info1 = F[10].split(';')
                info2 = F[11].split(';')
                for i in range(len(readids)):
                    tinfo1 = info1[i].split(',') 
                    tinfo2 = info2[i].split(',')
                    start1, end1, start2, end2 = int(tinfo1[0]), int(tinfo1[2]), int(tinfo2[0]), int(tinfo2[2])
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    if start1 <= start2:
                        readid2alignment[readids[i]].append((key, max(end1 - alignment_margin, 1), start2 + alignment_margin, '+', '*'))
                    else:
                        readid2alignment[readids[i]].append((key, max(end2 - alignment_margin, 1), start1 + alignment_margin, '-', '*'))

            elif mode == "d":
                size = F[10].split(';')
                info = F[11].split(';')
                for i in range(len(readids)):
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 1), min(tpos + alignment_margin, tlen - 1), tstrand, size[i]))

            elif mode == "i":
                size = F[10].split(';')
                info = F[11].split(';')
                median_size = round(statistics.median([int(x) for x in size if x not in ['-', '+']]))
                for i in range(len(readids)):
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    if size[i] in ['-', '+']: # breakpoint
                        if (tstrand == '+' and size[i] == '+') or (tstrand == '-' and size[i] == '-'):
                            readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 1), min(tpos + median_size + alignment_margin, tlen - 1), tstrand, size[i]))
                        else:
                            readid2alignment[readids[i]].append((key, max(tpos - median_size - alignment_margin, 1), min(tpos + alignment_margin, tlen - 1), tstrand, size[i]))
                    else: # Insertion identified by CIGAR
                        if tstrand == '+':
                            readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 1), min(tpos + int(size[i]) + alignment_margin, tlen - 1), tstrand, size[i]))
                        else:
                            readid2alignment[readids[i]].append((key, max(tpos - int(size[i]) - alignment_margin, 1), min(tpos + alignment_margin, tlen - 1), tstrand, size[i]))

    return(readid2alignment)

 
def gather_support_read_seq(rearrangement_file, insertion_file, deletion_file, output_file, tumor_bam, alignment_margin = 300):

    bamfile = pysam.AlignmentFile(tumor_bam, "rb")

    readid2alignment = {}
    set_readid2alignment(readid2alignment, rearrangement_file, 'r', alignment_margin)
    set_readid2alignment(readid2alignment, insertion_file, 'i', alignment_margin)
    set_readid2alignment(readid2alignment, deletion_file, 'd', alignment_margin)

    hout = open(output_file + ".tmp.unsorted", 'w')
    for read in bamfile.fetch():

        if read.query_name in readid2alignment and not read.is_secondary and not read.is_supplementary:

            for talignment in readid2alignment[read.query_name] :

                tkey, tstart, tend, tstrand, tsize = talignment

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(tstart - 1):tend]
                if tstrand == '-': part_seq = reverse_complement(part_seq)

                print('\t'.join([tkey, read.query_name, tsize, part_seq]), file = hout)

    hout.close()

    hout = open(output_file, 'w')
    subprocess.check_call(["sort", "-k1,1", output_file + ".tmp.unsorted"], stdout = hout)
    hout.close()
    os.remove(output_file + ".tmp.unsorted")



