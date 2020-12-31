#! /usr/bin/env python

import sys, os, subprocess, shutil 
from collections import Counter
import pysam
import parasail

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


def generate_paf_file(query_fasta, target_fasta, output_file):

    user_matrix = parasail.matrix_create("ACGT", 2, -2)

    with open(target_fasta, 'r') as hin:
        for line in hin:
            if line.startswith('>'): 
                tid = line.rstrip('\n').split(' ')[0].lstrip('>')
            else:
                tseq = line.rstrip('\n')

    with open(query_fasta, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            if line.startswith('>'):
                qid = line.rstrip('\n').lstrip('>')
            else:
                qseq = line.rstrip('\n')
                
                res = parasail.ssw(qseq, tseq, 3, 1, user_matrix)
                print("%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t*\t*\t60" %
                    (qid, len(qseq), res.read_begin1, res.read_end1,
                    tid, len(tseq), res.ref_begin1, res.ref_end1), file = hout)


def generate_racon_consensus(temp_key, tmp_dir):

    with open(tmp_dir + '/' + temp_key + ".tmp.seg.first.fa", 'w') as hout3: 
        subprocess.check_call(["head", "-n", "2", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout= hout3)
                     
    generate_paf_file(tmp_dir + '/' + temp_key + ".supporting_read.fa",
        tmp_dir + '/' + temp_key + ".tmp.seg.first.fa",
        tmp_dir + '/' + temp_key + ".parasail.paf")
                     
    with open(tmp_dir + '/' + temp_key + ".racon1.fa", 'w') as hout3:
        subprocess.check_call(["racon", "-u", 
            tmp_dir + '/' + temp_key + ".supporting_read.fa",
            tmp_dir + '/' + temp_key + ".parasail.paf",
            tmp_dir + '/' + temp_key + ".tmp.seg.first.fa"],
            stdout = hout3, stderr = subprocess.DEVNULL)

    with open(tmp_dir + '/' + temp_key + ".racon1.fa", 'r') as hin3, \
        open(tmp_dir + "/" + temp_key + ".racon1.mod.fa", 'w') as hout3:
        tid = hin3.readline()
        print(">temp_consensus", file = hout3)
        tseq = hin3.readline().rstrip('\n')
        print(tseq, file = hout3)

    generate_paf_file(tmp_dir + '/' + temp_key + ".supporting_read.fa",
        tmp_dir + '/' + temp_key + ".racon1.mod.fa",
        tmp_dir + '/' + temp_key + ".parasail2.paf")

    with open(tmp_dir + '/' + temp_key + ".racon2.fa", 'w') as hout3:
        subprocess.check_call(["racon", "-u",
            tmp_dir + '/' + temp_key + ".supporting_read.fa",
            tmp_dir + '/' + temp_key + ".parasail2.paf",
            tmp_dir + '/' + temp_key + ".racon1.mod.fa"],
            stdout = hout3, stderr = subprocess.DEVNULL)


def get_refined_bp(contig, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, h_log, rd_margin = 20, i_margin = 500):

    start1 = max(1, start1)
    start2 = max(1, start2)

    if mode != "i":

        bstart1, bend1 = max(int(start1) - rd_margin, 1), int(end1) + rd_margin
        bstart2, bend2 = max(int(start2) - rd_margin, 1), int(end2) + rd_margin

        region1_seq = fasta_file_ins.fetch(chr1, bstart1 - 1, bend1)
        region2_seq = fasta_file_ins.fetch(chr2, bstart2 - 1, bend2)

        if dir1 == '-': region1_seq = reverse_complement(region1_seq)
        if dir2 == '+': region2_seq = reverse_complement(region2_seq)

        sret = smith_waterman.sw_jump(contig, region1_seq, region2_seq)
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
            print("Alignment inconsistent!!", file = sys.stderr)

        print(score, contig_align, region1_align, region2_align, file = h_log)
        print(contig_seq, file = h_log)
        print(region_seq, file = h_log)

        return(bp_pos1, bp_pos2, inseq)
    
    else:
    
        contig_start = contig[:min(i_margin, len(contig))]
        contig_end = contig[-min(i_margin, len(contig)):]

        region_seq = fasta_file_ins.fetch(chr1, max(0, start1 - 100), end2 + 100)
        sret = smith_waterman.sw_jump(region_seq, contig_start, contig_end)
        if sret is None: return(None)
        score, region_align, contig_start_align, contig_end_align, region_seq, contig_seq = sret

        bp_pos1 = start1 - 100 + region_align[1]
        bp_pos2 = start1 - 100 + region_align[2]

        inseq_start = contig_start_align[1]
        inseq_end = len(contig) - (len(contig_end) - contig_end_align[0] + 1)
        inseq = contig[inseq_start:(inseq_end + 1)]

        print(score, region_align, contig_start_align, contig_end_align, file = h_log)
        print(region_seq, file = h_log)
        print(contig_seq, file = h_log)

        return(bp_pos1, bp_pos2, inseq)

        
def get_readid2alignment(input_file, mode, alignment_margin):

    readid2alignment = {}
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
                        readid2alignment[readids[i]].append((key, max(end1 - alignment_margin, 1), start2 + alignment_margin, '+'))
                    else:
                        readid2alignment[readids[i]].append((key, max(end2 - alignment_margin, 1), start1 + alignment_margin, '-'))

            elif mode == "d":
                size = F[10].split(';')
                info = F[11].split(';')
                for i in range(len(readids)):
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 1), min(tpos + alignment_margin, tlen - 1), tstrand))

            elif mode == "i":
                size = F[10].split(';')
                info = F[11].split(';')
                for i in range(len(readids)):
                    tinfo = info[i].split(',')
                    tpos, tlen, tstrand = int(tinfo[1]), int(tinfo[3]), tinfo[4]
                    if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                    if tstrand == "+":
                        readid2alignment[readids[i]].append((key, max(tpos - alignment_margin, 1), min(tpos + int(size[i]) + alignment_margin, tlen - 1), tstrand))
                    else:
                        readid2alignment[readids[i]].append((key, max(tpos - int(size[i]) - alignment_margin, 1), min(tpos + alignment_margin, tlen - 1), tstrand))

    return(readid2alignment)

 
def identify(rearrangement_file, insertion_file, deletion_file, output_file, tumor_bam, reference_fasta, debug, alignment_margin = 300):

    bamfile = pysam.AlignmentFile(tumor_bam, "rb")

    readid2alignment = get_readid2alignment(rearrangement_file, 'r', alignment_margin)
    readid2alignment.update(get_readid2alignment(insertion_file, 'i', alignment_margin))
    readid2alignment.update(get_readid2alignment(deletion_file, 'd', alignment_margin))

    hout = open(output_file + ".tmp.supporting_read.unsorted", 'w')
    for read in bamfile.fetch():

        if read.query_name in readid2alignment and not read.is_secondary and not read.is_supplementary:

            for talignment in readid2alignment[read.query_name]:

                tkey, tstart, tend, tstrand = talignment

                # if read.query_name == "cdca0470-c734-40e6-b0e8-6bfe60c41439":
                #     import pdb; pdb.set_trace()

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(tstart - 1):tend]
                if tstrand == '-': part_seq = reverse_complement(part_seq)
                # if read.is_reverse: part_seq = reverse_complement(part_seq)

                print('\t'.join([tkey, read.query_name, part_seq]), file = hout)

    hout.close()

    hout = open(output_file + ".tmp.supporting_read.sorted", 'w')
    subprocess.check_call(["sort", "-k1,1", output_file + ".tmp.supporting_read.unsorted"], stdout = hout)
    hout.close()
    subprocess.check_call(["rm", "-rf", output_file + ".tmp.supporting_read.unsorted"])
    
    fasta_file_ins = pysam.FastaFile(reference_fasta)

    tmp_dir = output_file + ".tmp_alignment_dir"
    if os.path.exists(tmp_dir): shutil.rmtree(tmp_dir)
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

                    """
                    hout3 = open(tmp_dir + '/' + temp_key + ".mafft_result.fa", 'w')
                    subprocess.check_call(["mafft", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout = hout3, stderr = subprocess.DEVNULL)
                    hout3.close()
                
                    tconsensus = get_consensus_from_mafft_result(tmp_dir + '/' + temp_key + ".mafft_result.fa")
                    print(temp_key + '\n' + tconsensus, file = hout_log)
                    # print(temp_key + '\n' + tconsensus)
                    """ 
                    generate_racon_consensus(temp_key, tmp_dir)
 
                    with open(tmp_dir + "/" + temp_key + ".racon2.fa") as hin2:
                        header = hin2.readline()
                        tconsensus = hin2.readline().rstrip('\n')
                    print(temp_key + '\n' + tconsensus, file = hout_log)

                    chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode = temp_key.split(',')
                    start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)       
                    bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, hout_log)
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

            """
            hout3 = open(tmp_dir + '/' + temp_key + ".mafft_result.fa", 'w')
            subprocess.check_call(["mafft", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout = hout3, stderr = subprocess.DEVNULL)
            hout3.close()

            tconsensus = get_consensus_from_mafft_result(tmp_dir + '/' + temp_key + ".mafft_result.fa")
            print(temp_key + '\n' + tconsensus, file = hout_log)
            """

            generate_racon_consensus(temp_key, tmp_dir)
                 
            with open(tmp_dir + "/" + temp_key + ".racon2.fa") as hin2:
                header = hin2.readline()
                tconsensus = hin2.readline().rstrip('\n')
            print(temp_key + '\n' + tconsensus, file = hout_log)

            chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode = temp_key.split(',')     
            start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2)
            bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, hout_log)
            if bret is not None:
                bp_pos1, bp_pos2, inseq = bret
                print(bp_pos1, bp_pos2, inseq, file = hout_log)
                print('', file = hout_log)
                print('\t'.join([chr1, str(bp_pos1), dir1, chr2, str(bp_pos2), dir2, inseq]), file = hout)

    hout_log.close()
    hout.close()

    if not debug:
        subprocess.check_call(["rm", "-rf", output_file + ".tmp.supporting_read.sorted"])
        shutil.rmtree(tmp_dir)


if __name__ == "__main__":

    """
    import sys
    rearrangement_file = sys.argv[1]
    deletion_file = sys.argv[2]
    output_file = sys.argv[3]
    tumor_bam = sys.argv[4]
    reference_fasta = sys.argv[5]

    identify(rearrangement_file, deletion_file, output_file, tumor_bam, reference_fasta)

    key = "X,104681514,104681542,+,X,104681515,104681543,-,i"
    tconsensus = "ATATGATTTAATAACCAATTCACAGAAAAGTTTAAATGACTGATAATAGCTACATGTGTTTCCTCTAGTAGGGCCATATTCATTTGGGTAGGTTACCAAAATCTTTAAATTAAAATAACGCATTTTCAACAGTAGGTGGGATGCCAAGCAACAACCAAAGGATAAAGGAGTTGGTAGGGATAAGAGTCAGAACAGCCTCCATTACATTACTGAAATAGAAGATATATTTGCTTTAAAAATTTTAATTTGGAGGTTGTATAAGAATATCACTTAGTCTATATGGTACTATAGACAAATATATATAGAAATGTAGGCATTTAGTGCTATAAATTTCCCTCTACACACTGCTTTGAATGGGTCCCAGAGATTCTGGTATGTGGTGTCTTTGTTCTCGTTGGTTTCAAAGAACATCTTTATTTCTGCCTTCATTTCAATATGTACCCAGTAGTCATTCAGGAGCAGGTTGTTCAGTTTCCAGTGTAGTTGAGCAGCTTTGGTGAGATTCTTAATCTGAGTTCTAGTTTGATTGCACTGTGGTCTGAGAGATAGTTTGTTATAATTTCTGTTCTTTTACATTTGCTGAGGAGATTTACTTCAACTATGTGGTCAATTTTGGAATAGGTGTGGTGTGGTGCTGAAAAAAATGTATATTCTGTTGATTTGGGGGAGAGTTCTGTAGATGTCTATTAGGTCTGCTTGGTGAGGGTTCAATCCTGGGTATCCTTGACTGGATTAAGAAAATGTGGCTACACCATGGAATACTATGCAGCCATAAAATGATGAGTTCATATCCTTTAAGGGACAATGAATGAAACATCATTCTCAGTAAACTATCATAAAACAAAAACAAAGCATATTCTACTCATAGGTGGGAATTGAACAATGAGTCACATGGACACAGGAGGAATATCACCTCTGGGGACTGTGGTGGGGTCGGGGGGGAGGATAGCATTGGGAATATACCTAATTAGATGACACATAGTGGGTGCAGCGCCAGCATGCACATGTATACATGTAACTAACCTGCCAATGTGCATGTGTACTAAAACTTAGGTATAATAAAAAAAAAAAAAAATAATAATAATAAAAAAAAAAAAAAAAAAAAAAAAAAAATATAGTAACACTCCCAAATACTGATGAATGTGCTCAGACCAGGAGAAATCTTACATGGAAAATAAAATAAGTTCTTGAAATCTTAACTTTAATAAATCATTTAGTTGATTGATTTTACATCAAAATTCTGAAACTAAACTCATTGTAAAATATAAGCGATATATATATATATATACCATATTTTGTCACTTTCTGAGAGCTGAGCAACTATAGTTGCATATAAATCACTTTAAAACTTATGTAATAGTGAGAATTCATTTTAAATTATAGTATTTTTA"
    """

    key = "9,12910928,12910950,+,9,12910929,12910951,-,i"
    tconsensus = "AACACATGAAGTATTTGGATCAGCTGGGTACAGTGGCTCACGCCTCTAATCTCAACACTTTGGGAGGCTGAGGCGGGTGGATCACTGAGGTCAGGAGTTTGAGACCAGCCTGGCCAACATAGTGAAGCTGTTTCTACTAAAATACAAAAATTAGCCAGGCGTGGCGGGCACCTGTAATCCCAGCTACTCAGGAGGCTAAGGCAGGAGAATCACTTGAACCTGGGAGGCAGAGGTTGCAGTGAGCCGAGATCGTGTCACTACACACTTCAACTTGGGTGACAGACTGAGATATACCTAATGCTAGATGACACATTAGTGGGTGCAGCGCACCAGCATAGCATGTATACATATGTAACTAACCTGCACAATGTGCACATGTCCCTAAAACTTAAGTATAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTATTTGAATCAAGGTCTTTTGAAGTCTCAGGCATTGACCATCTTGGTCTGTGAAGTCCTCACTATGATAAATATGCCAATTATGCAGCCTTCAAGTAGATCCTCAGTTACATGTACTGAGGCAGGAAAGGACTAAAGCCAAAAGCTCAATTACAGGCACCACTATGGGAGTCAAGCCCCTTATTTGAGGCTGAAACAAAATATCAGAAACCCCCAAGGATAATCCAGACTTCGAAATTAGCAGGCGATCAACTGAAAAGTCAGACATTGTTCTATAT"

    chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode = key.split(',')
    start1, end1, start2, end2 = int(start1), int(end1), int(start2), int(end2) 
    hout_log = open("test.txt", 'w')
    fasta_file_ins = pysam.FastaFile("/home/ubuntu/environment/workspace/seq_data/reference/GRCh37.fa")

    bret = get_refined_bp(tconsensus, fasta_file_ins, chr1, start1, end1, dir1, chr2, start2, end2, dir2, mode, hout_log)
    print(bret)

