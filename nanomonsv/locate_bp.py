#! /usr/bin/env python3

import sys, os, subprocess, shutil, statistics 
import pysam

from . import smith_waterman
from .my_seq import reverse_complement
from .logger import get_logger

logger = get_logger(__name__)


def get_refined_bp(contig, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, h_log, 
                   sw_jump_params, rd_margin = 20, i_margin = 500):

    start1 = max(1, start1)
    start2 = max(1, start2)

    if mode != "i":

        ref_len1 = fasta_file_ins.get_reference_length(chr1)
        ref_len2 = fasta_file_ins.get_reference_length(chr2)

        bstart1, bend1 = max(int(start1) - rd_margin, 1), int(end1) + rd_margin
        bstart2, bend2 = max(int(start2) - rd_margin, 1), int(end2) + rd_margin

        if ref_len1 < bend1: bend1 = ref_len1
        if ref_len2 < bend2: bend2 = ref_len2    

        region1_seq = fasta_file_ins.fetch(chr1, bstart1 - 1, bend1).upper()
        region2_seq = fasta_file_ins.fetch(chr2, bstart2 - 1, bend2).upper()
        

        if dir1 == '-': region1_seq = reverse_complement(region1_seq)
        if dir2 == '+': region2_seq = reverse_complement(region2_seq)

        sret = smith_waterman.sw_jump(contig, region1_seq, region2_seq,
                                      match_score = sw_jump_params[0], mismatch_penalty = sw_jump_params[1],
                                      gap_cost = sw_jump_params[2], jump_cost = sw_jump_params[3])

        if sret is None: return(None)
        score, contig_align, region1_align, region2_align, contig_seq, region_seq = sret

        bp_pos1 = bstart1 + region1_align[1] - 1 if dir1 == '+' else bend1 - region1_align[1] + 1 
        bp_pos2 = bstart2 + region2_align[0] - 1 if dir2 == '-' else bend2 - region2_align[0] + 1

        if contig_align[2] - contig_align[1] == 1:
            inseq = '---'
        elif contig_align[2] - contig_align[1] > 1:
            inseq = contig[(contig_align[1]):(contig_align[2] - 1)]
            if dir1 == '-': inseq = reverse_complement(inseq)
        else:
            logger.warning("Alignment inconsistent!!")

        print(score, contig_align, region1_align, region2_align, file = h_log)
        print(contig_seq, file = h_log)
        print(region_seq, file = h_log)

        return(bp_pos1, bp_pos2, inseq)
    
    else:
    
        contig_start = contig[:min(i_margin, len(contig))]
        contig_end = contig[-min(i_margin, len(contig)):]

        region_seq = fasta_file_ins.fetch(chr1, max(0, start1 - 100), end2 + 100)
        sret = smith_waterman.sw_jump(region_seq, contig_start, contig_end,
                                      match_score = sw_jump_params[0], mismatch_penalty = sw_jump_params[1], 
                                      gap_cost = sw_jump_params[2], jump_cost = sw_jump_params[3])

        if sret is None: return(None)
        score, region_align, contig_start_align, contig_end_align, region_seq, contig_seq = sret

        bp_pos1 = max(0, start1 - 100) + region_align[1]
        bp_pos2 = max(0, start1 - 100) + region_align[2]

        inseq_start = contig_start_align[1]
        inseq_end = len(contig) - (len(contig_end) - contig_end_align[0] + 1)
        inseq = contig[inseq_start:(inseq_end + 1)]

        print(score, region_align, contig_start_align, contig_end_align, file = h_log)
        print(region_seq, file = h_log)
        print(contig_seq, file = h_log)

        return(bp_pos1, bp_pos2, inseq)

        
 
def locate_bp(consensus_file, output_file, reference_fasta, sw_jump_params, debug):

    fasta_file_ins = pysam.FastaFile(reference_fasta)

    with open(consensus_file, 'r') as hin, open(output_file, 'w') as hout, \
        open(output_file + ".locate_bp.log", 'w') as hout_log:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            temp_key, tconsensus = F[0], F[1]
            print(temp_key + '\n' + tconsensus, file = hout_log)
            chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, cid = temp_key.split(',')
            start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)       
            bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, hout_log, sw_jump_params)
            if bret is not None:  
                bp_pos1, bp_pos2, inseq = bret 
                print(bp_pos1, bp_pos2, inseq, file = hout_log)
                print('', file = hout_log)
                print(f"{chr1}\t{bp_pos1}\t{dir1}\t{chr2}\t{bp_pos2}\t{dir2}\t{inseq}\t{mode}_{cid}", file = hout)
                # print('\t'.join([chr1, str(bp_pos1), dir1, chr2, str(bp_pos2), dir2, inseq]), file = hout)

    if not debug: os.remove(output_file + ".locate_bp.log")
