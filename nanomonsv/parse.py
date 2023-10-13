#! /usr/bin/env python3

import sys, os, subprocess, itertools
import pysam

from .logger import get_logger as logger
from .utils import get_alignment_object


def parse_alignment_info(input_alignment_file, reference_fasta, deletion_output_file, insertion_output_file, rearrangement_output_file,
    breakpoint_output_file, min_ins_size = 20, min_del_size = 30, min_clipping_size_for_bp = 200):

    """Parse BAM file to obtain putative SV supporting reads and their associated info."""

    hout_d = open(deletion_output_file, 'w')
    hout_i = open(insertion_output_file, 'w')
    hout_r = open(rearrangement_output_file, 'w')
    hout_b = open(breakpoint_output_file, 'w')

    alignment_h = get_alignment_object(input_alignment_file, reference_fasta)

    for read in alignment_h.fetch():

        # if read.is_secondary: continue

        query_name = read.query_name
        query_strand = '-' if read.is_reverse else '+'
        query_length = read.infer_read_length()

        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end
        mapping_quality = str(read.mapping_quality)
        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary

        cigar_stats = read.get_cigar_stats()
        num_M = cigar_stats[0][0]
        num_I = cigar_stats[0][1]
        num_D = cigar_stats[0][2]

        cigartuples = read.cigartuples
        left_hard_clipping_size, right_hard_clipping_size = 0, 0
        left_soft_clipping_size, right_soft_clipping_size = 0, 0
        if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]
        if cigartuples[0][0] == 4: left_soft_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] == 4: right_soft_clipping_size = cigartuples[-1][1]

        if not is_supplementary:
            if query_strand == '+':
                query_start = read.query_alignment_start + 1
                query_end = read.query_alignment_end
            else:
                query_start = query_length - read.query_alignment_end + 1
                query_end = query_length - read.query_alignment_start
        else:
            if query_strand == '+':
                query_start = left_hard_clipping_size + left_soft_clipping_size + 1
                query_end = query_length - right_hard_clipping_size - right_soft_clipping_size
            else:
                query_start = right_hard_clipping_size + right_soft_clipping_size + 1
                query_end = query_length - left_hard_clipping_size - left_soft_clipping_size


        query_pos_cur = query_start - 1 if query_strand == '+' else query_end
        query_pos_check = query_pos_cur
        reference_pos_cur = reference_start - 1
        reference_pos_check = reference_start - 1

        for cigar in cigartuples:
            if cigar[0] == 0:
                query_pos_cur = query_pos_cur + cigar[1] if query_strand == '+' else query_pos_cur - cigar[1]
                reference_pos_cur = reference_pos_cur + cigar[1]
            elif cigar[0] == 1:

                if cigar[1] >= min_ins_size:

                    tinfo = f'{query_start},{query_pos_cur},{query_end},{query_length},{query_strand},' + \
                        f'{mapping_quality},{num_M},{num_I - cigar[1]},{num_D},{is_supplementary},{is_secondary}'

                    print(f'{reference_name}\t{reference_pos_cur}\t{reference_pos_cur + 1}\t{query_name}\t{cigar[1]}\t+\t{tinfo}',
                        file = hout_i)

                query_pos_cur = query_pos_cur + cigar[1] if query_strand == '+' else query_pos_cur - cigar[1]

            elif cigar[0] == 2:
                if cigar[1] >= min_del_size:

                    """
                    if query_strand == '+':
                        query_start2 = query_pos_check + 1
                        query_end2 = query_pos_cur
                        reference_start2 = reference_pos_check + 1
                        reference_end2 = reference_pos_cur
                    else:
                        query_start2 = query_pos_cur + 1
                        query_end2 = query_pos_check
                        reference_start2 = reference_pos_check + 1
                        reference_end2 = reference_pos_cur
                    """

                    tinfo = f'{query_start},{query_pos_cur},{query_end},{query_length},{query_strand},' + \
                        f'{mapping_quality},{num_M},{num_I},{num_D - cigar[1]},{is_supplementary},{is_secondary}'

                    print(f'{reference_name}\t{reference_pos_cur}\t{reference_pos_cur + cigar[1]}\t{query_name}\t{cigar[1]}\t+\t{tinfo}',
                        file = hout_d)

                    reference_pos_cur = reference_pos_cur + cigar[1]
                else:
                    reference_pos_cur = reference_pos_cur + cigar[1]

        if query_strand == '+' and query_end != query_pos_cur:
            logger().error("query end inconsistent!! %s: %d != %d" % (query_name, query_end, query_pos_cur))
            sys.exit(1)
        if query_strand == '-' and query_start != query_pos_cur + 1:
            logger().error("query end inconsistent!! %s: %d != %d" % (query_name, query_end, query_pos_cur))
            sys.exit(1)

        """
        if query_strand == '+':
            query_start2 = query_pos_check + 1
            query_end2 = query_pos_cur
            reference_start2 = reference_pos_check + 1
            reference_end2 = reference_pos_cur
        else:
            query_start2 = query_pos_cur + 1
            query_end2 = query_pos_check
            reference_start2 = reference_pos_check + 1
            reference_end2 = reference_pos_cur

        import pdb; pdb.set_trace()
        if query_start != query_start2:
            print("query start inconsistent!!")
            sys.exit(1)
        if query_end != query_end2:
            print("query end inconsistent!!")
            sys.exit(1)
        if reference_start != reference_start2:
            print("reference start inconsistent!!")
            sys.exit(1)
        if reference_end != reference_end2:
            print("reference end inconsistent!!")
            sys.exit(1)
        """

        print(f'{query_name}\t{query_start}\t{query_end}\t{query_length}\t{query_strand}\t{reference_name}\t' + \
            f'{reference_start}\t{reference_end}\t{mapping_quality}\t{num_M}\t{num_I}\t{num_D}\t' + \
            f'{is_supplementary}\t{is_secondary}', file = hout_r)

        # print('\t'.join([query_name, str(query_start), str(query_end), str(query_length), query_strand, \
        #                  reference_name, str(reference_start), str(reference_end), mapping_quality, \
        #                  str(num_M), str(num_I), str(num_D), str(is_supplementary), str(is_secondary)]), file = hout_r)


        if left_soft_clipping_size + left_hard_clipping_size >= min_clipping_size_for_bp:

            if query_strand == '+':
                tinfo = f'{query_start},{query_start},{query_end},{query_length},{query_strand},' + \
                    f'{mapping_quality},{num_M},{num_I},{num_D},{is_supplementary},{is_secondary}'
            else:
                tinfo = f'{query_start},{query_end},{query_end},{query_length},{query_strand},' + \
                    f'{mapping_quality},{num_M},{num_I},{num_D},{is_supplementary},{is_secondary}'

            print(f'{reference_name}\t{reference_start - 1}\t{reference_start}\t{query_name}\t*\t-\t{tinfo}',
                file = hout_b)

        if right_soft_clipping_size + right_hard_clipping_size >= min_clipping_size_for_bp:

            if query_strand == '+':
                tinfo = f'{query_start},{query_end},{query_end},{query_length},{query_strand},' + \
                    f'{mapping_quality},{num_M},{num_I},{num_D},{is_supplementary},{is_secondary}'
            else:
                tinfo = f'{query_start},{query_start},{query_end},{query_length},{query_strand},' + \
                    f'{mapping_quality},{num_M},{num_I},{num_D},{is_supplementary},{is_secondary}'

            print(f'{reference_name}\t{reference_end - 1}\t{reference_end}\t{query_name}\t*\t+\t{tinfo}',
                file = hout_b)


    hout_d.close()
    hout_i.close()
    hout_r.close()
    hout_b.close()
    alignment_h.close()


def extract_bedpe_junction(input_file, output_file, split_alignment_check_margin1 = 50, split_alignment_check_margin2 = 50, minimum_ambiguity = 20):


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

            # if abs(qpos2[1] - qpos1[2]) <= split_alignment_check_margin:
            if qpos2[1] - qpos1[2] <= split_alignment_check_margin1 and qpos1[2] - qpos2[1] <= split_alignment_check_margin2:
                bp_flag = True
                tchr1, tstart1, tend1, tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1, tis_2nd1 = query2target[qpos1]
                tchr2, tstart2, tend2, tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2, tis_2nd2 = query2target[qpos2]

                if qpos2[1] - qpos1[2] > 0:
                    outward_ambiguity, inward_ambiguity = max(qpos2[1] - qpos1[2], minimum_ambiguity), minimum_ambiguity
                else:
                    outward_ambiguity, inward_ambiguity = minimum_ambiguity, max(qpos1[2] - qpos2[1], minimum_ambiguity)


                bchr1, bchr2 = tchr1, tchr2
                if qpos1[4] == '+':
                    bstart1, bend1, bstrand1 = max(int(tend1) - inward_ambiguity, 0), int(tend1) + outward_ambiguity, '+'
                else:
                    bstart1, bend1, bstrand1 = max(int(tstart1) - outward_ambiguity, 0), int(tstart1) + inward_ambiguity, '-'

                if qpos2[4] == '+':
                    bstart2, bend2, bstrand2 = max(int(tstart2) - outward_ambiguity, 0), int(tstart2) + inward_ambiguity, '-'
                else:
                    bstart2, bend2, bstrand2 = max(int(tend2) - inward_ambiguity, 0), int(tend2) + outward_ambiguity, '+'


                bread_name = qpos1[0]

                binfo1 = ','.join([str(qpos1[1]), '*', str(qpos1[2]), str(qpos1[3]), qpos1[4], tmapQ1, tnumM1, tnumI1, tnumD1, tis_supp1, tis_2nd1])
                binfo2 = ','.join([str(qpos2[1]), '*', str(qpos2[2]), str(qpos2[3]), qpos2[4], tmapQ2, tnumM2, tnumI2, tnumD2, tis_supp2, tis_2nd2])

                if bchr1 > bchr2 or (bchr1 == bchr2 and bstart1 > bstart2):
                    bchr1, bstart1, bend1, bstrand1, binfo1, bchr2, bstart2, bend2, bstrand2, binfo2 = \
                        bchr2, bstart2, bend2, bstrand2, binfo2, bchr1, bstart1, bend1, bstrand1, binfo1

                print(f'{bchr1}\t{bstart1}\t{bend1}\t{bchr2}\t{bstart2}\t{bend2}\t{bread_name}\t0\t{bstrand1}\t{bstrand2}\t{binfo1}\t{binfo2}',
                    file = hout)


    hout = open(output_file, 'w')

    temp_read_name = ''
    query2target = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            # skip secondary alignment
            # if F[13] == "True": continue

            if F[0] != temp_read_name:
                if temp_read_name != '' and len(query2target) > 1: print_bedpe_junction(query2target, hout)

                temp_read_name = F[0]
                query2target = {}

            query2target[(F[0], int(F[1]), int(F[2]), int(F[3]), F[4])] = (F[5], int(F[6]), int(F[7]), F[8], F[9], F[10], F[11], F[12], F[13])

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
