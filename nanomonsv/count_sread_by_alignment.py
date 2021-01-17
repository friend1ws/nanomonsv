#! /usr/bin/env python3

import sys, os, subprocess, shutil
import pysam
import parasail

from .pyssw import *
from .my_seq import reverse_complement

# function for gathering sequence read for realignment validation
def gather_local_read_for_realignment(sv_file, bam_file, output_file):

    bam_ps = pysam.AlignmentFile(bam_file, "rb")

    rname2key = {}
    key2rname2mapq = {}
    with open(sv_file, 'r') as hin:
        for line in hin:
            if line.startswith("#") or line.startswith("Chr_1"): continue

            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            # if tinseq == "---": tinseq = ''
            key = f"{tchr1},{tpos1},{tdir1},{tchr2},{tpos2},{tdir2},{tinseq}"

            if key not in key2rname2mapq: key2rname2mapq[key] = {}

            for read in bam_ps.fetch(tchr1, max(tpos1 - 100, 0), tpos1 + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

                if read.qname not in key2rname2mapq[key]: key2rname2mapq[key][read.qname] = [None, None]
                key2rname2mapq[key][read.qname][0] = read.mapping_quality

            for read in bam_ps.fetch(tchr2, max(tpos2 - 100, 0), tpos2 + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

                if read.qname not in key2rname2mapq[key]: key2rname2mapq[key][read.qname] = [None, None]
                key2rname2mapq[key][read.qname][1] = read.mapping_quality


    # remove duplicated keys
    for rname in rname2key:
        keys = list(set(rname2key[rname]))
        rname2key[rname] = keys

    with open(output_file + ".tmp.unsorted", 'w') as hout:
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
                    mapq1, mapq2 = key2rname2mapq[key][read.qname]
                    print(f"{key}\t{read.qname}\t{mapq1}\t{mapq2}\t{read_seq}", file = hout)

    with open(output_file, 'w') as hout:
        subprocess.call(["sort", "-k1,1", output_file + ".tmp.unsorted"], stdout = hout)

    os.remove(output_file + ".tmp.unsorted")
    bam_ps.close()


def ssw_check(query, target, use_ssw_lib):

    if use_ssw_lib:
        return(ssw_check_ssw_lib(query, target))
    else:
        return(ssw_check_parasail(query, target))


def ssw_check_ssw_lib(target, query):

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

            alignment_info[sQId] = [resPrint[0], rstart, rend, qstart, qend, strand]
        ssw.init_destroy(qProfile)
        ssw.init_destroy(qRcProfile)
       
    # return(supporting_reads)
    return(alignment_info)


def ssw_check_parasail(query, target):

    user_matrix = parasail.matrix_create("ACGT", 2, -2)

    alignment_info = {}
    for sQId, sQSeq, sQQual in read(query):

        sQSeq_r = reverse_complement(sQSeq)

        for sTId, sTSeq, STQual in read(target):

            res = parasail.ssw(sQSeq, sTSeq, 3, 1, user_matrix)
            res_r = parasail.ssw(sQSeq_r, sTSeq, 3, 1, user_matrix)

            if res.score1 > res_r.score1:
                score = res.score1
                qstart, qend = res.read_begin1 + 1., res.read_end1 + 1
                tstart, tend = res.ref_begin1 + 1, res.ref_end1 + 1
                strand = '+'
            else:
                score = res_r.score1
                qstart, qend = len(sQSeq) - res_r.read_end1, len(sQSeq) - res_r.read_begin1
                tstart, tend = res_r.ref_begin1 + 1, res_r.ref_end1 + 1
                strand = '-'

            alignment_info[sTId] = [score, int(qstart), int(qend), int(tstart), int(tend), strand]

    return(alignment_info)


class Alignment_counter(object):

    def __init__(self, output_count_file, output_alignment_info_file, reference_fasta,
        var_read_min_mapq, use_ssw_lib, debug, validate_sequence_length = 200,  
        score_ratio_thres = 1.2, start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20):

        self.temp_key = None
        self.temp_key2 = None # here, the insertion sequence is converted to its length. this is to shorten file names.
        self.hout_count = open(output_count_file, 'w')
        self.hout_ainfo = open(output_alignment_info_file, 'w')
        self.reference_fasta_h = pysam.FastaFile(reference_fasta)

        self.reference_segment_1 = None
        self.reference_segment_2 = None
        self.variant_segment_1 = None
        self.variant_segment_2 = None
        self.readid2mapq = {}

        self.validate_sequence_length = validate_sequence_length
        self.score_ratio_thres = score_ratio_thres
        self.start_pos_thres = start_pos_thres
        self.end_pos_thres = end_pos_thres
        self.var_ref_margin_thres = var_ref_margin_thres
        self.var_read_min_mapq = var_read_min_mapq
        self.use_ssw_lib = use_ssw_lib
        self.debug = debug

        self.tmp_dir = output_count_file + ".tmp_dir"
        os.makedirs(self.tmp_dir, exist_ok = True)


    def __del__(self):

        self.hout_count.close()
        self.hout_ainfo.close()
        self.reference_fasta_h.close()
        if not self.debug:
            shutil.rmtree(self.tmp_dir)

    def initialize(self, key):
        self.temp_key = key
        # the insertion sequence is converted to its length. this is to shorten file names.
        tkeys = self.temp_key.split(',')
        tkeys[-1] = '' if tkeys[-1] == '---' else tkeys[-1]
        tkeys[-1] = str(len(tkeys[-1]))
        self.temp_key2 = ','.join(tkeys)
        self.readid2mapq = {}
        self.generate_segment_fasta_files()
        self.temp_long_read_seq_file_h = open(self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", 'w')


    def generate_segment_fasta_files(self):

        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = self.temp_key.split(',')
        tpos1, tpos2 = int(tpos1), int(tpos2)
        tinseq = '' if tinseq == "---" else tinseq

        # reference_local_seq
        self.reference_segment_1 = self.reference_fasta_h.fetch(tchr1, max(tpos1 - self.validate_sequence_length - 1, 0), tpos1 + self.validate_sequence_length - 1)
        self.reference_segment_2 = self.reference_fasta_h.fetch(tchr2, max(tpos2 - self.validate_sequence_length - 1, 0), tpos2 + self.validate_sequence_length - 1)

        # variant_seq
        variant_seq = ""
        if tdir1 == '+':
            tseq = self.reference_fasta_h.fetch(tchr1, max(tpos1 - self.validate_sequence_length - 1, 0), tpos1 - 1)
        else:
            tseq = self.reference_fasta_h.fetch(tchr1, tpos1 - 1, tpos1 + self.validate_sequence_length - 1)
            tseq = reverse_complement(tseq)

        if tdir1 == "+":
            variant_seq = tseq + tinseq
        else:
            variant_seq = tseq + reverse_complement(tinseq)

        if tdir2 == '-':
            tseq = self.reference_fasta_h.fetch(tchr2, tpos2 - 1, tpos2 + self.validate_sequence_length - 1)
        else:
            tseq = self.reference_fasta_h.fetch(tchr2, max(tpos2 - self.validate_sequence_length - 1, 0), tpos2 - 1)
            tseq = reverse_complement(tseq)

        variant_seq = variant_seq + tseq

        self.variant_segment_1 = variant_seq[:min(2 * self.validate_sequence_length, len(variant_seq))]
        self.variant_segment_2 = variant_seq[-min(2 * self.validate_sequence_length, len(variant_seq)):]


    def add_long_read_seq(self, treadid, tseq, tmapq1, tmapq2):

        self.readid2mapq[treadid] = [int(tmapq1) if tmapq1 != "None" else None, int(tmapq2) if tmapq2 != "None" else None]
        print(f">{treadid}\n{tseq}", file = self.temp_long_read_seq_file_h)


    def count_alignment_and_print(self):

        def is_short_del_dup(key):
            keys = key.split(',')
            if keys[6] == "---": keys[6] == ''
            if keys[0] == keys[3] and keys[2] == '+' and keys[5] == '-' and int(keys[4]) - int(keys[1]) + len(keys[6]) < 100:
                return(True)
            elif keys[0] == keys[3] and keys[2] == '-' and keys[5] == '+' and int(keys[4]) - int(keys[1]) + len(keys[6]) < 100:
                return(True)
            else:
                return(False)

        # close file hundle for writing long_read_seq.fa if it is open
        if self.temp_long_read_seq_file_h is not None: self.temp_long_read_seq_file_h.close()

        with open(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.variant_segment_1, file = hout)
        with open(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_2.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.variant_segment_2, file = hout)
        with open(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_1.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.reference_segment_1, file = hout)
        with open(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_2.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.reference_segment_2, file = hout)

        alignment_info_var_1 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)
        alignment_info_var_2 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_2.fa", 
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib) 

        # if is_short_del_dup(self.temp_key):
        alignment_info_ref_1 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_1.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)
        alignment_info_ref_2 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_2.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)

        all_rnames = list(set(list(alignment_info_var_1.keys()) + list(alignment_info_var_2.keys())))

        supporting_reads = [rname for rname in all_rnames if \
            (alignment_info_var_1[rname][0] > self.score_ratio_thres * len(self.variant_segment_1) and \
            alignment_info_var_1[rname][1] < self.start_pos_thres * len(self.variant_segment_1) and \
            alignment_info_var_1[rname][2] > self.end_pos_thres * len(self.variant_segment_1)) or \
            (alignment_info_var_2[rname][0] > self.score_ratio_thres * len(self.variant_segment_2) and \
            alignment_info_var_2[rname][1] < self.start_pos_thres * len(self.variant_segment_2) and \
            alignment_info_var_2[rname][2] > self.end_pos_thres * len(self.variant_segment_2))]

        # for short deletion or insertion
        # if is_short_del_dup(self.temp_key):
        supporting_reads = [rname for rname in supporting_reads if \
            (alignment_info_var_1[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres and \
            alignment_info_var_1[rname][0] >= alignment_info_ref_2[rname][0] + self.var_ref_margin_thres) or \
            (alignment_info_var_2[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres and \
             alignment_info_var_2[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres)]

        # filtering by mapping qualities
        supporting_reads = [rname for rname in supporting_reads if \
            self.readid2mapq[rname][0] is not None and self.readid2mapq[rname][0] >= self.var_read_min_mapq and \
            self.readid2mapq[rname][1] is not None and self.readid2mapq[rname][1] >= self.var_read_min_mapq]


        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = self.temp_key.split(',')
        print(f"{tchr1}\t{tpos1}\t{tdir1}\t{tchr2}\t{tpos2}\t{tdir2}\t{tinseq}\t{len(all_rnames)}\t{len(supporting_reads)}", file = self.hout_count)
        
        sread_info = { rname: alignment_info_var_1[rname] + alignment_info_var_2[rname] for rname in supporting_reads}
        for rname in supporting_reads:
            sinfo = '\t'.join([str(x) for x in alignment_info_var_1[rname] + alignment_info_var_2[rname]])
            print(f"{tchr1}\t{tpos1}\t{tdir1}\t{tchr2}\t{tpos2}\t{tdir2}\t{tinseq}\t{rname}\t{sinfo}", file = self.hout_ainfo)


def count_sread_by_alignment(sv_file, bam_file, output_count_file, output_alignment_info_file, reference_fasta, 
    var_read_min_mapq, use_ssw_lib, debug):

    gather_local_read_for_realignment(sv_file, bam_file, output_count_file + ".tmp.local_read_for_realignment")
 
    alignment_counter = Alignment_counter(output_count_file, output_alignment_info_file, reference_fasta,
        var_read_min_mapq, use_ssw_lib, debug)

    with open(output_count_file + ".tmp.local_read_for_realignment", 'r') as hin:
        for line in hin:
            tkey, treadid, tmapq1, tmpaq2, tseq = line.rstrip('\n').split('\t')
            if alignment_counter.temp_key != tkey:
                if alignment_counter.temp_key is not None:
                    alignment_counter.count_alignment_and_print()
                alignment_counter.initialize(tkey)

            alignment_counter.add_long_read_seq(treadid, tseq, tmapq1, tmpaq2)

        if alignment_counter.temp_key is not None:
            alignment_counter.count_alignment_and_print()


    if not debug:
        os.remove(output_count_file + ".tmp.local_read_for_realignment")


