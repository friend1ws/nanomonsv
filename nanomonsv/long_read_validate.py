#! /usr/bin/env python

import os, tempfile, subprocess, shutil
import pysam
from onebreak.long_read_validate import ssw_check
from onebreak.my_seq import reverse_complement

def long_read_validate_by_alignment(sv_file, output_file, bam_file, reference, debug, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8, var_ref_margin_thres = 10):

    bam_ps = pysam.AlignmentFile(bam_file, "rb")

    rname2key = {}
    with open(sv_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')               
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2])

            for read in bam_ps.fetch(tchr1, max(tpos1 - 100, 0), tpos1 + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

            for read in bam_ps.fetch(tchr2, max(tpos2 - 100, 0), tpos2 + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

    # remove duplicated keys
    for rname in rname2key:
        keys = list(set(rname2key[rname]))
        rname2key[rname] = keys

    hout = open(output_file + ".tmp3.long_read_seq.unsorted", 'w')
    for read in bam_ps.fetch():

        flags = format(int(read.flag), "#014b")[:1:-1]
    
        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads             
        if flags[10] == "1": continue

        if read.qname in rname2key:

            read_seq = read.query_sequence
            if flags[4] == "1": read_seq = reverse_complement(read_seq)

            for key in rname2key[read.qname]:

                print(key + '\t' + read.qname + '\t' + read_seq, file = hout)

    hout.close()
    bam_ps.close()

    hout = open(output_file + ".tmp3.long_read_seq.sorted", 'w')
    subprocess.call(["sort", "-k1,1", output_file + ".tmp3.long_read_seq.unsorted"], stdout = hout)
    hout.close()
    subprocess.call(["rm" ,"-rf", output_file + ".tmp3.long_read_seq.unsorted"])

    # my_seq.get_seq function could be used. But this procedure is repeatead many times and using pysam class may be good for the IO.
    reference_fasta = pysam.FastaFile(os.path.abspath(reference))
 
    key2contig = {}
    with open(sv_file, 'r') as hin:

        for line in hin:
            F = line.rstrip('\n').split('\t')

            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2])
            if tinseq == "---": tinseq = ''
 
            # reference_local_seq
            # reference_local_seq_1 = reference_fasta.fetch(tchr1, max(tpos1 - validate_sequence_length - 1, 0), tpos1 + validate_sequence_length - 1) 
            # reference_local_seq_2 = reference_fasta.fetch(tchr2, max(tpos2 - validate_sequence_length - 1, 0), tpos2 + validate_sequence_length - 1)  

            # variant_seq
            variant_seq = ""
            if tdir1 == '+':
                tseq = reference_fasta.fetch(tchr1, max(tpos1 - validate_sequence_length - 1, 0), tpos1 - 1)
            else:  
                tseq = reference_fasta.fetch(tchr1, tpos1 - 1, tpos1 + validate_sequence_length - 1)
                tseq = reverse_complement(tseq)

            variant_seq = tseq + tinseq
 
            if tdir2 == '-':
                tseq = reference_fasta.fetch(tchr2, tpos2 - 1, tpos2 + validate_sequence_length - 1)
            else:
                tseq = reference_fasta.fetch(tchr2, max(tpos2 - validate_sequence_length - 1, 0), tpos2 - 1)
                tseq = reverse_complement(tseq)

            variant_seq = variant_seq + tseq

            variant_seq_1 = variant_seq[:min(2 * validate_sequence_length, len(variant_seq))]
            variant_seq_2 = variant_seq[-min(2 * validate_sequence_length, len(variant_seq)):]

            # print(key)
            # key2contig[key] = [variant_seq, reference_local_seq_1,reference_local_seq_2]
            key2contig[key] = [variant_seq_1, variant_seq_2]

       
    # import pdb; pdb.set_trace()

    tmp_dir = tempfile.mkdtemp()
    # tmp_dir = "tmp"
    # os.makedirs(tmp_dir)
    # print(tmp_dir)

    temp_key = ""
    temp_id2seq = {}
    temp_junc_seq = ""
    temp_total_read_count = 0
    key2sread = {}
    key2sread_info = {}
    key2sread_count = {}
    key2sread_count_all = {}
    with open(output_file + ".tmp3.long_read_seq.sorted") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if temp_key != F[0]:
                if temp_key != "":
                    hout2.close()
                    alignment_info_var_1 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_1.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
                    alignment_info_var_2 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_2.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
                    all_keys = list(set(alignment_info_var_1))

                    supporting_read_keys = [key for key in all_keys if \
                        (alignment_info_var_1[key][0] > score_ratio_thres * len(variant_seq_1) and \
                        alignment_info_var_1[key][1] < start_pos_thres * len(variant_seq_1) and \
                        alignment_info_var_1[key][2] > end_pos_thres * len(variant_seq_1)) or \
                        (alignment_info_var_2[key][0] > score_ratio_thres * len(variant_seq_2) and \
                        alignment_info_var_2[key][1] < start_pos_thres * len(variant_seq_2) and \
                        alignment_info_var_2[key][2] > end_pos_thres * len(variant_seq_2))]

                    key2sread[temp_key] = supporting_read_keys
                    key2sread_info[temp_key] = { key: alignment_info_var_1[key] + alignment_info_var_2[key] for key in supporting_read_keys}
                    key2sread_count[temp_key] = len(supporting_read_keys)
                    key2sread_count_all[temp_key] = len(all_keys)

                    # print(temp_key + '\t' + str(key2sread_count_all[temp_key]) + '\t' + str(key2sread_count[temp_key]))

                hout2 = open(tmp_dir + '/' + F[0] + ".long_read_seq.fa", 'w')
                temp_key = F[0]
                temp_total_read_count = 0
                FF = temp_key.split(',')
                variant_seq_1, variant_seq_2 = key2contig[temp_key]
                with open(tmp_dir + '/' + F[0] + ".variant_seq_1.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + variant_seq_1, file = hout1)
                with open(tmp_dir + '/' + F[0] + ".variant_seq_2.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + variant_seq_2, file = hout1)

            print('>' + F[1] + '\n' + F[2], file = hout2)
            temp_total_read_count = temp_total_read_count + 1


        # last transaction
        if temp_key != "":
            hout2.close()

            alignment_info_var_1 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_1.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
            alignment_info_var_2 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_2.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
            all_keys = list(set(alignment_info_var_1))

            supporting_read_keys = [key for key in all_keys if \
                (alignment_info_var_1[key][0] > score_ratio_thres * len(variant_seq_1) and \
                alignment_info_var_1[key][1] < start_pos_thres * len(variant_seq_1) and \
                alignment_info_var_1[key][2] > end_pos_thres * len(variant_seq_1)) or \
                (alignment_info_var_2[key][0] > score_ratio_thres * len(variant_seq_2) and \
                alignment_info_var_2[key][1] < start_pos_thres * len(variant_seq_2) and \
                alignment_info_var_2[key][2] > end_pos_thres * len(variant_seq_2))]

            key2sread[temp_key] = supporting_read_keys
            key2sread_info[temp_key] = { key: alignment_info_var_1[key] + alignment_info_var_2[key] for key in supporting_read_keys}
            key2sread_count[temp_key] = len(supporting_read_keys)
            key2sread_count_all[temp_key] = len(all_keys)


    shutil.rmtree(tmp_dir)
    if not debug: subprocess.call(["rm" ,"-rf", output_file + ".tmp3.long_read_seq.sorted"])

    return([key2sread_count, key2sread_count_all, key2sread_info])



def validate_main(result_file, tumor_bam, output, sread_file, reference, control_bam, debug):

    key2sread_count_tumor, key2sread_count_all_tumor, key2sread_info_tumor = long_read_validate_by_alignment(result_file, output, tumor_bam, reference, debug, score_ratio_thres = 1.3, start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20)

    if control_bam is not None:
        key2sread_count_control, key2sread_count_all_control, key2sread_info_control = long_read_validate_by_alignment(result_file, output, control_bam, reference, debug, score_ratio_thres = 1.3, start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20)

    hout = open(output, 'w')
    with open(result_file, 'r') as hin:

        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2])

            sread_count_tumor = key2sread_count_tumor[key] if key in key2sread_count_tumor else 0
            sread_count_all_tumor = key2sread_count_all_tumor[key] if key in key2sread_count_all_tumor else 0
            # sread_id = ','.join(key2sread_tumor[key])

            if control_bam is not None:
                sread_count_control = key2sread_count_control[key] if key in key2sread_count_control else 0
                sread_count_all_control = key2sread_count_all_control[key] if key in key2sread_count_all_control else 0

            if control_bam is not None:
                print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor) + '\t' + \
                      str(sread_count_all_control) + '\t' + str(sread_count_control), file = hout)
                      # str(sread_count_all_control) + '\t' + str(sread_count_control) + '\t' + sread_id, file = hout)
            else:
                # print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor) + '\t' + sread_id, file = hout)
                print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor), file = hout)

    hout.close()

    hout = open(sread_file, 'w')
    for key in key2sread_info_tumor:
        for sread in key2sread_info_tumor[key]:
            sinfo = '\t'.join([str(x) for x in key2sread_info_tumor[key][sread]])
            print(key + '\t' + sread + '\t' + sinfo, file = hout)

    hout.close()

if __name__ == "__main__":

    import sys
    result_file = sys.argv[1]
    tumor_bam = sys.argv[2]
    output = sys.argv[3]
    reference = sys.argv[4]
    control_bam = sys.argv[5]

    validate_main(result_file, tumor_bam, output, reference, control_bam)


