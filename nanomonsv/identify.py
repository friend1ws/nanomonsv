#! /usr/bin/env python

import sys, os, subprocess
from collections import Counter
import pysam

from . import smith_waterman

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))



def get_consensus_from_mafft_result(input_file):

    id2seq = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                tid = line
                id2seq[tid] = ''
            else:
                id2seq[tid] = id2seq[tid] + line


    ind2bases = {}
    for tid in id2seq:
        seq = id2seq[tid]
        for i in range(len(seq)):
            if i not in ind2bases: ind2bases[i] = []
            ind2bases[i].append(seq[i])


    seq_len = len(list(ind2bases))
    consensus = ''
    for i in range(seq_len):
        # import pdb; pdb.set_trace()
        mycounter = Counter(ind2bases[i] )
        consensus = consensus + mycounter.most_common()[0][0]

    consensus = consensus.replace('-', '').upper()

    return(consensus)



def get_refined_bp(contig, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, h_log):

    start1 = max(1, start1)
    start2 = max(1, start2)

    region1_seq = fasta_file_ins.fetch(chr1, start1 - 1, end1)
    region2_seq = fasta_file_ins.fetch(chr2, start2 - 1, end2)

    if dir1 == '-': region1_seq = reverse_complement(region1_seq)
    if dir2 == '+': region2_seq = reverse_complement(region2_seq)

    sret = smith_waterman.sw_jump(contig, region1_seq, region2_seq)
    if sret is None: return(None)
    score, contig_align, region1_align, region2_align, contig_seq, region_seq = sret

    bp_pos1 = start1 + region1_align[1] - 1 if dir1 == '+' else end1 - region1_align[1] + 1 
    bp_pos2 = start2 + region2_align[0] - 1 if dir2 == '-' else end2 - region2_align[0] + 1

    if contig_align[2] - contig_align[1] == 1:
        inseq = '---'
    elif contig_align[2] - contig_align[1] > 1:
        inseq = contig[(contig_align[1]):(contig_align[2] - 1)]
    else:
        print("Alignment consistent!!", file = sys.stderr)

    print(score, contig_align, region1_align, region2_align, file = h_log)
    print(contig_seq, file = h_log)
    print(region_seq, file = h_log)

    return(bp_pos1, bp_pos2, inseq)


def get_readid2alignment(input_file, mode, alignment_margin):

    readid2alignment = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = ','.join([F[0], F[1], F[2], F[8], F[3], F[4], F[5], F[9]])
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
                        readid2alignment[readids[i]].append((key, max(end1 - alignment_margin, 0), start2 + alignment_margin, '+'))
                    else:
                        readid2alignment[readids[i]].append((key, max(end2 - alignment_margin, 0), start1 + alignment_margin, '-'))

            elif mode == "d":
                size = F[10].split(';')
                info = F[11].split(';')
                for i in range(len(readids)):
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 0), min(tpos + alignment_margin, tlen - 1), tstrand))

    return(readid2alignment)

 
def identify(rearrangement_file, deletion_file, output_file, tumor_bam, reference_fasta, alignment_margin = 300):

    bamfile = pysam.AlignmentFile(tumor_bam, "rb")

    """
    readid2alignment = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            key = ','.join([F[0], F[1], F[2], F[8], F[3], F[4], F[5], F[9]]) 

            readids = F[6].split(';')

            if mode == "r":
                info1 = F[10].split(';')
                info2 = F[11].split(';')
            elif mode == "d":
                size = F[10].split(';')
                info = F[11].split(';')

            for i in range(len(readids)):

                if mode == "r": 
                    tinfo1 = info1[i].split(',')
                    tinfo2 = info2[i].split(',')
                    start1, end1, start2, end2 = int(tinfo1[0]), int(tinfo1[2]), int(tinfo2[0]), int(tinfo2[2])
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    if start1 <= start2:
                        readid2alignment[readids[i]].append((key, end1 - alignment_margin, start2 + alignment_margin, '+'))
                    else:
                        readid2alignment[readids[i]].append((key, end2 - alignment_margin, start1 + alignment_margin, '-'))

                elif mode == "d":
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 0), min(tpos + alignment_margin, tlen - 1), tstrand))
    """

    readid2alignment = get_readid2alignment(rearrangement_file, 'r', alignment_margin)
    readid2alignment.update(get_readid2alignment(deletion_file, 'd', alignment_margin))

    hout = open(output_file + ".tmp.supporting_read.unsorted", 'w')
    for read in bamfile.fetch():

        if read.query_name in readid2alignment and not read.is_secondary and not read.is_supplementary:

            for talignment in readid2alignment[read.query_name]:

                tkey, tstart, tend, tstrand = talignment

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(tstart - 1):tend]
                if tstrand == '-': part_seq = reverse_complement(part_seq)
                # if read.is_reverse: part_seq = reverse_complement(part_seq)

                print('\t'.join([tkey, read.query_name, part_seq]), file = hout)

    hout.close()

    hout = open(output_file + ".tmp.supporting_read.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp.supporting_read.unsorted"], stdout = hout)
    hout.close()

    fasta_file_ins = pysam.FastaFile(reference_fasta)

    tmp_dir = "tmp"
    os.makedirs(tmp_dir)
 
    temp_key = ''
    # hout = open(output_file + ".tmp.consensus.fastq", 'w') 
    hout_log = open(tmp_dir + "/consensus_alignment.log", 'w')

    hout = open(output_file, 'w') 
    with open(output_file + ".tmp.supporting_read.sorted", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if temp_key != F[0]: 

                if temp_key != '':
                    hout2.close()

                    hout3 = open(tmp_dir + '/' + temp_key + ".mafft_result.fa", 'w')
                    subprocess.check_call(["mafft", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout = hout3, stderr = subprocess.DEVNULL)
                    hout3.close()

                    tconsensus = get_consensus_from_mafft_result(tmp_dir + '/' + temp_key + ".mafft_result.fa")
                    print(temp_key + '\n' + tconsensus, file = hout_log)
                    print(temp_key)
 
                    chr1, start1, end1, dir1, chr2, start2, end2, dir2 = temp_key.split(',')
                    start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)       
                    bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, hout_log)
                    if bret is not None:  
                        bp_pos1, bp_pos2, inseq = bret 
                        print(bp_pos1, bp_pos2, inseq, file = hout_log)
                        print('', file = hout_log)
                        print('\t'.join([chr1, str(bp_pos1), dir1, chr2, str(bp_pos2), dir2, inseq]), file = hout)
 
                temp_key = F[0]          
                hout2 = open(tmp_dir + '/' + temp_key + ".supporting_read.fa", 'w')

            print('>' + F[1] + '\n' + F[2], file = hout2)

        # last treatment
        if temp_key != '':
            hout2.close()

            hout3 = open(tmp_dir + '/' + temp_key + ".mafft_result.fa", 'w')
            subprocess.check_call(["mafft", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout = hout3, stderr = subprocess.DEVNULL)
            hout3.close()

            tconsensus = get_consensus_from_mafft_result(tmp_dir + '/' + temp_key + ".mafft_result.fa")
            print(temp_key + '\n' + tconsensus, file = hout_log)

            chr1, start1, end1, dir1, chr2, start2, end2, dir2 = temp_key.split(',')     
            start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)
            bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, hout_log)
            if bret is not None:
                bp_pos1, bp_pos2, inseq = bret
                print(bp_pos1, bp_pos2, inseq, file = hout_log)
                print('', file = hout_log)
                print('\t'.join([chr1, str(bp_pos1), dir1, chr2, str(bp_pos2), dir2, inseq]), file = hout)

    hout_log.close()
    hout.close()
            


if __name__ == "__main__":

    import sys
    rearrangement_file = sys.argv[1]
    deletion_file = sys.argv[2]
    output_file = sys.argv[3]
    tumor_bam = sys.argv[4]
    reference_fasta = sys.argv[5]

    identify(rearrangement_file, deletion_file, output_file, tumor_bam, reference_fasta)

 
