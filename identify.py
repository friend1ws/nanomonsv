#! /usr/bin/env python

import os, subprocess
from collections import Counter
import pysam

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



def identify(input_file, output_file, tumor_bam, alignment_margin = 300):

    """
    bamfile = pysam.AlignmentFile(tumor_bam, "rb")

    readid2alignment = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            key = ','.join([F[0], F[1], F[2], F[8], F[3], F[4], F[5], F[9]]) 

            readids = F[6].split(';')
            info1 = F[10].split(';')
            info2 = F[11].split(';')

            for i in range(len(readids)):
                tinfo1 = info1[i].split(',')
                tinfo2 = info2[i].split(',')
                start1, end1, start2, end2 = int(tinfo1[0]), int(tinfo1[1]), int(tinfo2[0]), int(tinfo2[1])
                if readids[i] not in readid2alignment: readid2alignment[readids[i]] = []
                if start1 <= start2:
                    readid2alignment[readids[i]].append((key, end1 - alignment_margin, start2 + alignment_margin, '+'))
                else:
                    readid2alignment[readids[i]].append((key, end2 - alignment_margin, start1 + alignment_margin, '-'))


    hout = open(output_file + "tmp.supporting_read.unsorted", 'w')
    for read in bamfile.fetch():

        if read.query_name in readid2alignment and not read.is_secondary and not read.is_supplementary:

            for talignment in readid2alignment[read.query_name]:

                tkey, tstart, tend, tstrand = talignment

                read_seq = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence
                part_seq = read_seq[(tstart - 1):tend]
                if tstrand == '-': part_seq = reverse_complement(part_seq)
                # if read.is_reverse: part_seq = reverse_complement(part_seq)

                print('\t'.join([tkey, read.query_name, part_seq]), file = hout1)

    hout.close()
    """

    hout = open(output_file + ".tmp.supporting_read.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp.supporting_read.unsorted"], stdout = hout)
    hout.close()

    
    tmp_dir = "tmp"
    os.makedirs(tmp_dir)
 
    temp_key = ''
    hout = open(output_file + ".tmp.consensus.fastq", 'w') 
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
                    print('@' + temp_key, file = hout)
                    print(tconsensus, file = hout)
                    print('+', file = hout)
                    print(''.join(['I' for x in range(len(tconsensus))]), file = hout)
               
                    print(temp_key)
 
                temp_key = F[0]          
                hout2 = open(tmp_dir + '/' + temp_key + ".supporting_read.fa", 'w')

            print('>' + F[1] + '\n' + F[2], file = hout2)

        # last treatment
        if temp_key != '':
            hout2.close()

            hout3 = open(tmp_dir + '/' + temp_key + ".mafft_result.fa", 'w')
            subprocess.check_call(["mafft", tmp_dir + '/' + temp_key + ".supporting_read.fa"], stdout = hout3)
            hout3.close()

            tconsensus = get_consensus_from_mafft_result(tmp_dir + '/' + temp_key + ".mafft_result.fa")
            print('@' + temp_key, file = hout)
            print(tconsensus, file = hout)
            print('+', file = hout)
            print(''.join(['I' for x in range(len(contig))]), file = hout)

    hout.close()
            


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    tumor_bam = sys.argv[3]

    identify(input_file, output_file, tumor_bam)
 
