#! /usr/bin/env python

import subprocess, itertools
import pysam

def parse_alignment_info(input_bam, output_file):

    hout = open(output_file, 'w') 
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch():

        # print(read)
        if read.is_secondary: continue

        query_name = read.query_name
        query_strand = '-' if read.is_reverse else '+'
        query_length = read.infer_read_length()

        reference_name = read.reference_name
        reference_start = str(read.reference_start + 1)
        reference_end = str(read.reference_end)
        mapping_quality = str(read.mapping_quality)
        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary
        
        cigar_stats = read.get_cigar_stats()
        num_M = cigar_stats[0][0]
        num_I = cigar_stats[0][1]
        num_D = cigar_stats[0][2]

        cigartuples = read.cigartuples
        left_hard_clipping_size, right_hard_clipping_size = 0, 0
        if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]

        if not is_supplementary:
            if query_strand == '+':
                query_start = read.query_alignment_start + 1 
                query_end = read.query_alignment_end
            else:
                query_start = query_length - read.query_alignment_end + 1
                query_end = query_length - read.query_alignment_start
        else:
            if query_strand == '+':
                query_start = left_hard_clipping_size + 1
                query_end = query_length - right_hard_clipping_size
            else:
                query_start = right_hard_clipping_size + 1
                query_end = query_length - left_hard_clipping_size
            
        print('\t'.join([query_name, str(query_start), str(query_end), str(query_length), query_strand, \
                         reference_name, reference_start, reference_end, mapping_quality, \
                         str(num_M), str(num_I), str(num_D), str(is_supplementary), str(is_secondary)]), file = hout)

    hout.close()
    bamfile.close()


def extract_bedpe_junction(input_file, output_file, check_margin = 50, minimum_ambiguity = 20):

   
    def print_bedpe_junction(query2target, hout):

        query_list = list(query2target)
        for qpos_comb in list(itertools.combinations(query_list, 2)):

            if qpos_comb[0][1] < qpos_comb[1][1]:
                qpos1, qpos2 = qpos_comb[0], qpos_comb[1]
            else:
                qpos1, qpos2 = qpos_comb[1], qpos_comb[0]
            
            # if the first region completely covers the second region
            if qpos1[2] >= qpos2[2]: continue
            
            # if there is a significant overlap         
            if (qpos1[2] - qpos2[1]) / (qpos2[2] - qpos1[2]) >= 0.2: continue
            
            if abs(qpos2[1] - qpos1[2]) <= check_margin:
                bp_flag = True
                tchr1, tstart1, tend1, tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1 = query2target[qpos1]
                tchr2, tstart2, tend2, tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2 = query2target[qpos2]
                
                if qpos2[1] - qpos1[2] > 0:
                    outward_ambiguity, inward_ambiguity = minimum_ambiguity, max(qpos2[1] - qpos1[2], minimum_ambiguity)
                else:
                    outward_ambiguity, inward_ambiguity = max(qpos1[2] - qpos2[1], minimum_ambiguity), minimum_ambiguity
                
                bchr1, bchr2 = tchr1, tchr2
                if qpos1[4] == '+': 
                    bstart1, bend1, bstrand1 = int(tend1) - inward_ambiguity, int(tend1) + outward_ambiguity, '+'
                else:
                    bstart1, bend1, bstrand1 = int(tstart1) - outward_ambiguity, int(tstart1) + inward_ambiguity, '-'
                
                if qpos2[4] == '+': 
                    bstart2, bend2, bstrand2 = int(tstart2) - outward_ambiguity, int(tstart2) + inward_ambiguity, '-'
                else:
                    bstart2, bend2, bstrand2 = int(tend2) - inward_ambiguity, int(tend2) + outward_ambiguity, '+'
                
                bread_name = qpos1[0]
                
                binfo1 = ','.join([str(qpos1[1]), str(qpos1[2]), tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1])
                binfo2 = ','.join([str(qpos2[1]), str(qpos2[2]), tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2])
                
                if bchr1 > bchr2 or (bchr1 == bchr2 and bstart1 > bstart2):
                    bchr1, bstart1, bend1, bstrand1, binfo1, bchr2, bstart2, bend2, bstrand2, binfo2 = \
                        bchr2, bstart2, bend2, bstrand2, binfo2, bchr1, bstart1, bend1, bstrand1, binfo1
                
                print('\t'.join([bchr1, str(bstart1), str(bend1), bchr2, str(bstart2), str(bend2), bread_name, "0",
                                 bstrand1, bstrand2, binfo1, binfo2]), file = hout)

 
    hout = open(output_file, 'w')

    temp_read_name = ''
    query2target = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[11] == "True": continue

            if F[0] != temp_read_name:
                if temp_read_name != '' and len(query2target) > 1: print_bedpe_junction(query2target, hout)
                            
                temp_read_name = F[0]
                query2target = {}

            query2target[(F[0], int(F[1]), int(F[2]), int(F[3]), F[4])] = (F[5], int(F[6]), int(F[7]), F[8], F[9], F[10], F[11], F[12])

    if temp_read_name != '' and len(query2target) > 1: print_bedpe_junction(query2target, hout)

    hout.close()

 
if __name__ == "__main__":

    import sys

    input_bam = sys.argv[1]
    output_prefix = sys.argv[2]

    parse_alignment_info(input_bam, output_prefix + ".tmp.alignment_info.txt")

    hout = open(output_prefix + ".tmp.alignment_info.name_sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", output_prefix + ".tmp.alignment_info.txt"], stdout = hout)
    hout.close()

    extract_bedpe_junction(output_prefix + ".tmp.alignment_info.name_sorted.txt", output_prefix + ".tmp.junction.bedpe")

    hout = open(output_prefix + ".tmp.junction.sorted.bedpe", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", "-k5,5n", "-k6,6n", output_prefix + ".tmp.junction.bedpe"], stdout = hout)
    hout.close()
  
    hout = open(output_prefix + ".junction.sorted.bedpe.gz", 'w')
    subprocess.check_call(["bgzip", "-f", "-c", output_prefix + ".tmp.junction.sorted.bedpe"], stdout = hout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", output_prefix + ".junction.sorted.bedpe.gz"])




