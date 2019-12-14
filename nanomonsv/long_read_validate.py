#! /usr/bin/env python

import os, tempfile, subprocess, shutil
import pysam

from .pyssw import *
from .my_seq import reverse_complement
# from onebreak.long_read_validate import ssw_check
# from onebreak.my_seq import reverse_complement

def ssw_check(target, query):

    nMatch = 2
    nMismatch = 2
    nOpen = 3
    nExt = 1
    sLibPath = ""
 
    lEle = []
    dRc = {} 
    dEle2Int = {}
    dInt2Ele = {}
    # init DNA score matrix
    lEle = ['A', 'C', 'G', 'T', 'N']
    dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N': 'N', 'a':'T', 'c':'G', 'g':'C', 't':'A', 'n': 'N'} 
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = nMatch
            else:
                lScore[i*nEleNum+j] = -nMismatch

    # translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore

    # set flag
    nFlag = 1

    # check whether libssw.so is in LD_LIBRARY_PATH
    sLibPath = ""
    for ld_path in os.environ["LD_LIBRARY_PATH"].split(':'):
        # print ld_path
        if os.path.exists(ld_path + "/libssw.so"):
            sLibPath = ld_path # + "/libssw.so"
            break
    if sLibPath == "":
        print("cannot find libssw.so in LD_LIBRARY_PATH", file = sys.stderr)
        sys.exit(1)

    ssw = ssw_lib.CSsw(sLibPath)
    # supporting_reads = []
    alignment_info = {}

    # iterate query sequence
    for sQId,sQSeq,sQQual in read(query):

        sQSeq = sQSeq.strip('\n')

        # build query profile
        qNum = to_int(sQSeq, lEle, dEle2Int)
        qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

        # for reverse complement
        # import pdb; pdb.set_trace()
        sQRcSeq = ''.join([dRc[x] for x in sQSeq[::-1]])
        qRcNum = to_int(sQRcSeq, lEle, dEle2Int)
        qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)

        # set mask len
        if len(sQSeq) > 30:
            nMaskLen = int(len(sQSeq) / 2)
        else:
            nMaskLen = 15

        # iter target sequence
        for sRId,sRSeq,_ in read(target):

            rNum = to_int(sRSeq, lEle, dEle2Int)

            # format ofres: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
            res = align_one(ssw, qProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)

            # align rc query
            resRc = align_one(ssw, qRcProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    
            # build cigar and trace back path
            if resRc == None or res[0] > resRc[0]:
                resPrint = res
                rstart, rend = resPrint[2] + 1, resPrint[3] + 1
                qstart, qend = resPrint[4] + 1, resPrint[5] + 1
                strand = '+'
                sCigar, sQ, sA, sR = buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
            else:
                resPrint = resRc
                rstart, rend = resPrint[2] + 1, resPrint[3] + 1
                qstart, qend = len(sQSeq) - resPrint[5], len(sQSeq) - resPrint[4] 
                strand = '-'
                sCigar, sQ, sA, sR = buildPath(sQRcSeq, sRSeq, resRc[4], resRc[2], resRc[8])

            # import pdb; pdb.set_trace()
            # if int(resPrint[0]) > score_ratio_thres * len(sRSeq) and int(resPrint[2]) + 1 < start_pos_thres * len(sRSeq) and int(resPrint[3]) + 1 > end_pos_thres * len(sRSeq):
            # supporting_reads.append([sQId, resPrint[0], resPrint[2] + 1, resPrint[3] + 1])
            # alignment_info[sQId] = [resPrint[0], resPrint[2] + 1, resPrint[3] + 1, resPrint[4] + 1, resPrint[5] + 1, strand]
            alignment_info[sQId] = [resPrint[0], rstart, rend, qstart, qend, strand]
        ssw.init_destroy(qProfile)
        ssw.init_destroy(qRcProfile)
       
    # return(supporting_reads)
    return(alignment_info)


def long_read_validate_by_alignment(sv_file, output_file, bam_file, reference, debug, validate_sequence_length = 200, score_ratio_thres = 1.4, start_pos_thres = 0.2, end_pos_thres = 0.8, var_ref_margin_thres = 10):
    
    def is_short_del_dup(key):
        keys = key.split(',')
        if keys[6] == "---": keys[6] == ''
        if keys[0] == keys[3] and keys[2] == '+' and keys[5] == '-' and int(keys[4]) - int(keys[1]) + len(keys[6]) < 100:
            return(True)
        elif keys[0] == keys[3] and keys[2] == '-' and keys[5] == '+' and int(keys[4]) - int(keys[1]) + len(keys[6]) < 100:
            return(True)
        else:
            return(False)

    
    bam_ps = pysam.AlignmentFile(bam_file, "rb")

    rname2key = {}
    with open(sv_file, 'r') as hin:
        for line in hin:
            if line.startswith("#") or line.startswith("Chr_1"): continue 

            F = line.rstrip('\n').split('\t')               
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            if tinseq == "---": tinseq = ''
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2, str(len(tinseq))])

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
            if line.startswith("#") or line.startswith("Chr_1"): continue

            F = line.rstrip('\n').split('\t')

            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            if tinseq == "---": tinseq = ''
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2, str(len(tinseq))])
 
            # reference_local_seq
            reference_local_seq_1 = reference_fasta.fetch(tchr1, max(tpos1 - validate_sequence_length - 1, 0), tpos1 + validate_sequence_length - 1) 
            reference_local_seq_2 = reference_fasta.fetch(tchr2, max(tpos2 - validate_sequence_length - 1, 0), tpos2 + validate_sequence_length - 1)  

            # variant_seq
            variant_seq = ""
            if tdir1 == '+':
                tseq = reference_fasta.fetch(tchr1, max(tpos1 - validate_sequence_length - 1, 0), tpos1 - 1)
            else:  
                tseq = reference_fasta.fetch(tchr1, tpos1 - 1, tpos1 + validate_sequence_length - 1)
                tseq = reverse_complement(tseq)

            if tdir1 == "+":
                variant_seq = tseq + tinseq
            else:
                variant_seq = tseq + reverse_complement(tinseq)

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
            # if add_reference_alignment:
            #     key2contig[key] = [variant_seq_1, variant_seq_2, reference_local_seq_1,reference_local_seq_2]
            # else:
            #     key2contig[key] = [variant_seq_1, variant_seq_2]
            key2contig[key] = [variant_seq_1, variant_seq_2, reference_local_seq_1,reference_local_seq_2]

       
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
                    if is_short_del_dup(temp_key):
                        alignment_info_ref_1 = ssw_check(tmp_dir + '/' + temp_key + ".reference_local_seq_1.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
                        alignment_info_ref_2 = ssw_check(tmp_dir + '/' + temp_key + ".reference_local_seq_2.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")

                    all_keys = list(set(list(alignment_info_var_1.keys()) + list(alignment_info_var_2.keys())))

                    if is_short_del_dup(temp_key):
                        supporting_read_keys = [key for key in all_keys if \
                            (alignment_info_var_1[key][0] > score_ratio_thres * len(variant_seq_1) and \
                            alignment_info_var_1[key][1] < start_pos_thres * len(variant_seq_1) and \
                            alignment_info_var_1[key][2] > end_pos_thres * len(variant_seq_1) and \
                            alignment_info_var_1[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres and \
                            alignment_info_var_1[key][0] >= alignment_info_ref_2[key][0] + var_ref_margin_thres) or \
                            (alignment_info_var_2[key][0] > score_ratio_thres * len(variant_seq_2) and \
                            alignment_info_var_2[key][1] < start_pos_thres * len(variant_seq_2) and \
                            alignment_info_var_2[key][2] > end_pos_thres * len(variant_seq_2) and \
                            alignment_info_var_2[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres and \
                            alignment_info_var_2[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres)]
                    else:
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
    
                """
                if add_reference_alignment:
                    variant_seq_1, variant_seq_2, reference_local_seq_1, reference_local_seq_2 = key2contig[temp_key]
                    with open(tmp_dir + '/' + F[0] + ".variant_seq_1.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + variant_seq_1, file = hout1) 
                    with open(tmp_dir + '/' + F[0] + ".variant_seq_2.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + variant_seq_2, file = hout1)
                    with open(tmp_dir + '/' + F[0] + ".reference_local_seq_1.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + reference_local_seq_1, file = hout1)
                    with open(tmp_dir + '/' + F[0] + ".reference_local_seq_2.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + reference_local_seq_2, file = hout1)

                else:
                    variant_seq_1, variant_seq_2 = key2contig[temp_key]
                    with open(tmp_dir + '/' + F[0] + ".variant_seq_1.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + variant_seq_1, file = hout1)
                    with open(tmp_dir + '/' + F[0] + ".variant_seq_2.fa", 'w') as hout1:
                        print('>' + F[0] + '\n' + variant_seq_2, file = hout1)
                """

                variant_seq_1, variant_seq_2, reference_local_seq_1, reference_local_seq_2 = key2contig[temp_key]
                with open(tmp_dir + '/' + F[0] + ".variant_seq_1.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + variant_seq_1, file = hout1) 
                with open(tmp_dir + '/' + F[0] + ".variant_seq_2.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + variant_seq_2, file = hout1)
                with open(tmp_dir + '/' + F[0] + ".reference_local_seq_1.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + reference_local_seq_1, file = hout1) 
                with open(tmp_dir + '/' + F[0] + ".reference_local_seq_2.fa", 'w') as hout1:
                    print('>' + F[0] + '\n' + reference_local_seq_2, file = hout1)

            print('>' + F[1] + '\n' + F[2], file = hout2)
            temp_total_read_count = temp_total_read_count + 1


        # last transaction
        if temp_key != "":
            hout2.close()

            """
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
            """

            alignment_info_var_1 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_1.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
            alignment_info_var_2 = ssw_check(tmp_dir + '/' + temp_key + ".variant_seq_2.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
            if is_short_del_dup(temp_key):
                alignment_info_ref_1 = ssw_check(tmp_dir + '/' + temp_key + ".reference_local_seq_1.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")
                alignment_info_ref_2 = ssw_check(tmp_dir + '/' + temp_key + ".reference_local_seq_2.fa", tmp_dir + '/' + temp_key + ".long_read_seq.fa")

            all_keys = list(set(alignment_info_var_1))

            if is_short_del_dup(temp_key):
                supporting_read_keys = [key for key in all_keys if \
                    (alignment_info_var_1[key][0] > score_ratio_thres * len(variant_seq_1) and \
                    alignment_info_var_1[key][1] < start_pos_thres * len(variant_seq_1) and \
                    alignment_info_var_1[key][2] > end_pos_thres * len(variant_seq_1) and \
                    alignment_info_var_1[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres and \
                    alignment_info_var_1[key][0] >= alignment_info_ref_2[key][0] + var_ref_margin_thres) or \
                    (alignment_info_var_2[key][0] > score_ratio_thres * len(variant_seq_2) and \
                    alignment_info_var_2[key][1] < start_pos_thres * len(variant_seq_2) and \
                    alignment_info_var_2[key][2] > end_pos_thres * len(variant_seq_2) and \
                    alignment_info_var_2[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres and \
                    alignment_info_var_2[key][0] >= alignment_info_ref_1[key][0] + var_ref_margin_thres)]
            else:   
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



def long_read_validate_main(result_file, tumor_bam, output, sread_file, reference, control_bam, debug):

    key2sread_count_tumor, key2sread_count_all_tumor, key2sread_info_tumor = long_read_validate_by_alignment(result_file, output, tumor_bam, reference, debug, score_ratio_thres = 1.2, start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20)

    if control_bam is not None:
        key2sread_count_control, key2sread_count_all_control, key2sread_info_control = long_read_validate_by_alignment(result_file, output, control_bam, reference, debug, score_ratio_thres = 1.2, start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20)

    hout = open(output, 'w')
    with open(result_file, 'r') as hin:

        for line in hin:
            if line.startswith("#") or line.startswith("Chr_1"): continue

            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            if tinseq == "---": tinseq = ''
            key = ','.join([tchr1, str(tpos1), tdir1, tchr2, str(tpos2), tdir2, str(len(tinseq))])

            sread_count_tumor = key2sread_count_tumor[key] if key in key2sread_count_tumor else 0
            sread_count_all_tumor = key2sread_count_all_tumor[key] if key in key2sread_count_all_tumor else 0
            # sread_id = ','.join(key2sread_tumor[key])

            if control_bam is not None:
                sread_count_control = key2sread_count_control[key] if key in key2sread_count_control else 0
                sread_count_all_control = key2sread_count_all_control[key] if key in key2sread_count_all_control else 0

            if control_bam is not None:
                print('\t'.join(F[:7]) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor) + '\t' + \
                      str(sread_count_all_control) + '\t' + str(sread_count_control), file = hout)
                      # str(sread_count_all_control) + '\t' + str(sread_count_control) + '\t' + sread_id, file = hout)
            else:
                # print('\t'.join(F) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor) + '\t' + sread_id, file = hout)
                print('\t'.join(F[:7]) + '\t' + str(sread_count_all_tumor) + '\t' + str(sread_count_tumor), file = hout)

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

    long_read_validate_main(result_file, tumor_bam, output, output + ".sread.txt", reference, control_bam, False)



