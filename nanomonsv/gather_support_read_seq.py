#! /usr/bin/env python3

import sys, os, subprocess, statistics
import pysam

from .logger import get_logger as logger
from .my_seq import reverse_complement
from .utils import get_alignment_object


def set_readid2alignment(readid2alignment, input_file, mode, alignment_margin):

    # readid2alignment = {}
    cid = 0
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = f"{F[0]},{F[1]},{F[2]},{F[8]},{F[3]},{F[4]},{F[5]},{F[9]},{mode},{cid}"
                # ','.join([F[0], F[1], F[2], F[8], F[3], F[4], F[5], F[9], mode])
            # readids = F[6].split(';')
            readids = eval(F[6])

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

            cid = cid + 1



def set_readid2alignment_sbnd(readid2alignment_sbnd, input_file, alignment_margin):

    cid = 0
    with open(input_file, 'r') as hin:
        for line in hin:
            tchr, tstart, tend, treadids, _, tstrand, tinfos = line.rstrip('\n').split('\t')
            tkey = f"{tchr},{tstart},{tend},{tstrand},b,{cid}" # ','.join([tchr, tstart, tend, tstrand])

            # readids = treadids.split(';')
            readids = eval(treadids)

            infos = tinfos.split(';')
            for i in range(len(readids)):
                ttinfo = infos[i].split(',')
                qpos, qlen, qstrand = int(ttinfo[1]), int(ttinfo[3]), ttinfo[4]
                if readids[i] not in readid2alignment_sbnd: readid2alignment_sbnd[readids[i]] = []
                if tstrand == qstrand:
                    readid2alignment_sbnd[readids[i]].append((tkey, max(1, qpos - alignment_margin), qlen, qstrand, '*'))
                else:
                    readid2alignment_sbnd[readids[i]].append((tkey, 1, min(qpos + alignment_margin, qlen), qstrand, '*'))

            cid = cid + 1


def gather_support_read_seq(rearrangement_file, insertion_file, deletion_file, output_file, tumor_alignment_file, reference_fasta,
    single_breakend_file = None, output_file_sbind = None, alignment_margin = 300):

    is_sbnd = True if single_breakend_file is not None else False


    alignment_h = get_alignment_object(tumor_alignment_file, reference_fasta)

    readid2alignment = {}
    set_readid2alignment(readid2alignment, rearrangement_file, 'r', alignment_margin)
    set_readid2alignment(readid2alignment, insertion_file, 'i', alignment_margin)
    set_readid2alignment(readid2alignment, deletion_file, 'd', alignment_margin)

    hout = open(output_file + ".tmp.unsorted", 'w')

    if is_sbnd:
        readid2alignment_sbnd = {}
        set_readid2alignment_sbnd(readid2alignment_sbnd, single_breakend_file, alignment_margin)
        hout_sbnd = open(output_file + ".sbnd.tmp.unsorted", 'w')

    for read in alignment_h.fetch():

        if read.is_secondary or read.is_supplementary: continue

        if read.query_name in readid2alignment:

            for talignment in readid2alignment[read.query_name]:

                akey, astart, aend, astrand, asize = talignment

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(astart - 1):aend]
                if astrand == '-': part_seq = reverse_complement(part_seq)

                print('\t'.join([akey, read.query_name, asize, part_seq]), file = hout)


        if is_sbnd and read.query_name in readid2alignment_sbnd:

            for talignment in readid2alignment_sbnd[read.query_name]:

                akey, astart, aend, astrand, asize = talignment
                tchr, tstart, tend, tstrand, _, tcid = akey.split(',')

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(astart - 1):aend]
                if astrand != tstrand: part_seq = reverse_complement(part_seq)

                print('\t'.join([akey, read.query_name, asize, part_seq]), file = hout_sbnd)


    alignment_h.close()
    hout.close()
    if is_sbnd: hout_sbnd.close()

    with open(output_file, 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", output_file + ".tmp.unsorted"], stdout = hout)
    os.remove(output_file + ".tmp.unsorted")

    if is_sbnd:
        with open(output_file_sbind, 'w') as hout:
            subprocess.check_call(["sort", "-k1,1", output_file + ".sbnd.tmp.unsorted"], stdout = hout)
        os.remove(output_file + ".sbnd.tmp.unsorted")
