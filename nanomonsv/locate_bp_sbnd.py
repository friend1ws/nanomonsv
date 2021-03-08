#! /usr/bin/env python3

import sys, os, subprocess, shutil, statistics 
import pysam
import parasail

from .my_seq import reverse_complement

from .logger import get_logger

logger = get_logger(__name__)

       
def get_refined_bp_sbnd(tconsensus, fasta_file_h, tchr, tstart, tend, tdir, hout_log, margin = 200):

    tconsensus_part = tconsensus[:1000] if len(tconsensus) > 1000 else tconsensus

    ref_len = fasta_file_h.get_reference_length(tchr)
    if tstart < 1: tstart = 1
    if ref_len < tend: tend = ref_len 

    if tdir == '+':
        qseq = fasta_file_h.fetch(tchr, max(int(tstart) - margin, 0), int(tend))
    else:
        qseq = fasta_file_h.fetch(tchr, max(int(tstart) - 1, 0), int(tend) + margin)
        qseq = reverse_complement(qseq)

    user_matrix = parasail.matrix_create("ACGT", 1, -2)
    res = parasail.ssw(qseq, tconsensus, 3, 1, user_matrix)
    if res is None:
        logger.debug(f"Alignment for breakpoint localization failed for {tchr},{tstart},{tend},{tdir}")
        return None


    if tdir == '+':
        bp_pos_reference = tend - (len(qseq) - res.read_end1 - 1)
    else:
        bp_pos_reference = tstart + (len(qseq) - res.read_end1 - 1) 

    tconsensus_after = tconsensus[(res.ref_end1 + 1):]

    return (bp_pos_reference, tconsensus_after)

 
def locate_bp_sbnd(consensus_file, output_file, reference_fasta, debug):

    fasta_file_ins = pysam.FastaFile(reference_fasta)

    with open(consensus_file, 'r') as hin, open(output_file, 'w') as hout, \
        open(output_file + ".locate_bp.sbnd.log", 'w') as hout_log:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            temp_key, tconsensus = F[0], F[1]
            print(temp_key + '\n' + tconsensus, file = hout_log)
            tchr, tstart, tend, tdir, _, tcid = temp_key.split(',')
            tstart, tend = int(tstart), int(tend)  
            bret = get_refined_bp_sbnd(tconsensus, fasta_file_ins, tchr, tstart, tend, tdir, hout_log)
            if bret is not None:  
                bp_pos, tconsensus_after = bret
                print(bp_pos, tconsensus_after, file = hout_log)
                print('', file = hout_log)
                print(f"{tchr}\t{bp_pos}\t{tdir}\t{tconsensus_after}\tb_{tcid}", file = hout)
                # print('\t'.join([tchr, str(bp_pos), tdir, tconsensus_after]), file = hout)

    if not debug: os.remove(output_file + ".locate_bp.sbnd.log")


