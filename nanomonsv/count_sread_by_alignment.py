#! /usr/bin/env python3

import sys, os, subprocess, shutil, random, copy
import pysam
import parasail

from .pyssw import *
from .my_seq import reverse_complement
from .utils import get_alignment_object

# when the total read number exceeds max_read_num, then max_read_num reads are randomly selected
def bam_subsample_fetch(alignment_h, tchr, tstart, tend, max_read_num = 500):

    rec = 0
    for read in alignment_h.fetch(tchr, tstart, tend):
        rec = rec + 1

    selected_inds = random.sample(range(rec), min(max_read_num, rec))

    rec2 = 0
    for read in alignment_h.fetch(tchr, tstart, tend):
        if rec2 in selected_inds:
            yield read
        rec2 = rec2 + 1

# function for gathering sequence read for realignment validation
def gather_local_read_for_realignment(sv_file, alignment_file, output_file, reference_fasta,
    sbnd_file = None, output_file_sbnd = None, validate_sequence_length = 200, check_read_max_num = 500,
    sort_option = None):

    alignment_h = get_alignment_object(alignment_file, reference_fasta)

    rname2key = {}
    key2rname2mapq = {}
    with open(sv_file, 'r') as hin:
        for line in hin:
            if line.startswith("#") or line.startswith("Chr_1"): continue
            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6], F[7]
            # if tinseq == "---": tinseq = ''
            key = f"{tchr1},{tpos1},{tdir1},{tchr2},{tpos2},{tdir2},{tinseq},{tid}"

            if key not in key2rname2mapq: key2rname2mapq[key] = {}

            for read in bam_subsample_fetch(alignment_h, tchr1, max(tpos1 - 100, 0), tpos1 + 100):
    
                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

                if read.qname not in key2rname2mapq[key]: key2rname2mapq[key][read.qname] = [None, None]
                key2rname2mapq[key][read.qname][0] = read.mapping_quality

            for read in bam_subsample_fetch(alignment_h, tchr2, max(tpos2 - 100, 0), tpos2 + 100):

                if read.qname not in rname2key: rname2key[read.qname] = []
                rname2key[read.qname].append(key)

                if read.qname not in key2rname2mapq[key]: key2rname2mapq[key][read.qname] = [None, None]
                key2rname2mapq[key][read.qname][1] = read.mapping_quality

    # remove duplicated keys
    for rname in rname2key:
        keys = list(set(rname2key[rname]))
        rname2key[rname] = keys

    # for sbnd candidate
    if sbnd_file is not None:
        rname2key_sbnd = {}
        key2rname2mapq_sbnd = {}
        with open(sbnd_file, 'r') as hin:
            for line in hin:
                if line.startswith("#") or line.startswith("Chr_1"): continue
                F = line.rstrip('\n').split('\t')
                tchr, tpos, tdir, tseq, tid = F[0], int(F[1]), F[2], F[3][:validate_sequence_length], F[4]
                key = f"{tchr},{tpos},{tdir},{tseq},{tid}"

                if key not in key2rname2mapq_sbnd: key2rname2mapq_sbnd[key] = {}

                for read in bam_subsample_fetch(alignment_h, tchr, max(tpos - 100, 0), tpos + 100):

                    if read.qname not in rname2key_sbnd: rname2key_sbnd[read.qname] = []
                    rname2key_sbnd[read.qname].append(key)
                    
                    key2rname2mapq_sbnd[key][read.qname] = read.mapping_quality

        # remove duplicated keys
        for rname in rname2key_sbnd:
            keys = list(set(rname2key_sbnd[rname]))
            rname2key_sbnd[rname] = keys

 
    hout = open(output_file + ".tmp.unsorted", 'w')
    if sbnd_file is not None: hout_sbnd = open(output_file_sbnd + ".tmp.unsorted", 'w')
    for read in alignment_h.fetch():

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

        if sbnd_file is not None and read.qname in rname2key_sbnd:

            read_seq = read.query_sequence
            if flags[4] == "1": read_seq = reverse_complement(read_seq)

            for key in rname2key_sbnd[read.qname]:
                mapq = key2rname2mapq_sbnd[key][read.qname]
                print(f"{key}\t{read.qname}\t{mapq}\tNone\t{read_seq}", file = hout_sbnd)

    hout.close()
    if sbnd_file is not None: hout_sbnd.close()

    with open(output_file, 'w') as hout:
        subprocess.call(["sort", "-k1,1"] + sort_option.split(" ") + [output_file + ".tmp.unsorted"], stdout = hout)
    os.remove(output_file + ".tmp.unsorted")

    if sbnd_file is not None:
        with open(output_file_sbnd, 'w') as hout:
            subprocess.call(["sort", "-k1,1"] + sort_option.split(" ") + [output_file_sbnd + ".tmp.unsorted"], stdout = hout)
        os.remove(output_file_sbnd + ".tmp.unsorted")

    alignment_h.close()


def ssw_check(query, target, use_ssw = False):

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
        var_read_min_mapq, score_ratio_thres, use_ssw_lib, debug, validate_sequence_length = 200,  
        start_pos_thres = 0.1, end_pos_thres = 0.9, var_ref_margin_thres = 20):

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

    def initialize(self, key, is_sbnd = False):
        self.temp_key = key
        # the insertion sequence is converted to its length. this is to shorten file names.
        tkeys = self.temp_key.split(',')
        tkeys[-2] = '' if tkeys[-2] == '---' else tkeys[-1]
        tkeys[-2] = str(len(tkeys[-2]))
        self.temp_key2 = ','.join(tkeys)

        self.is_inseq = True if tkeys[-2] != '' else False

        if is_sbnd == False:
            self.is_short_del_dup = False
            if tkeys[6] == "---": tkeys[6] = ''
            if tkeys[0] == tkeys[3] and tkeys[2] == '+' and tkeys[5] == '-' and int(tkeys[4]) - int(tkeys[1]) + len(tkeys[6]) < 300:
                self.is_short_del_dup = True
            elif tkeys[0] == tkeys[3] and tkeys[2] == '-' and tkeys[5] == '+' and int(tkeys[4]) - int(tkeys[1]) + len(tkeys[6]) < 300:
                self.is_short_del_dup = True

        self.readid2mapq = {}
        self.temp_long_read_seq_file_h = open(self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", 'w')
        if is_sbnd:
            self.generate_segment_fasta_files_sbnd()
        else:
            self.generate_segment_fasta_files()

    def generate_segment_fasta_files(self):

        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid = self.temp_key.split(',')
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


    def generate_segment_fasta_files_sbnd(self):

        tchr, tpos, tdir, tcontig, tid = self.temp_key.split(',')
        tpos = int(tpos)

        # referene_local_seq
        self.reference_segment_1 = self.reference_fasta_h.fetch(tchr, max(tpos - self.validate_sequence_length - 1, 0), tpos + self.validate_sequence_length - 1)

        # variant_seq
        variant_seq = ''
        if tdir == '+':
            tseq = self.reference_fasta_h.fetch(tchr, max(tpos - self.validate_sequence_length - 1, 0), tpos)
            tseq = tseq + tcontig
        else:
            tseq = self.reference_fasta_h.fetch(tchr, tpos - 1, tpos + self.validate_sequence_length - 1)
            tseq = reverse_complement(tseq)
        self.variant_segment_1 = tseq + tcontig


    def add_long_read_seq(self, treadid, tseq, tmapq1, tmapq2):

        self.readid2mapq[treadid] = [int(tmapq1) if tmapq1 != "None" else None, int(tmapq2) if tmapq2 != "None" else None]
        print(f">{treadid}\n{tseq}", file = self.temp_long_read_seq_file_h)


    def count_alignment_and_print(self):

        # close file hundle for writing long_read_seq.fa if it is open
        if self.temp_long_read_seq_file_h is not None: self.temp_long_read_seq_file_h.close()

        with open(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.variant_segment_1, file = hout)
        with open(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_2.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.variant_segment_2, file = hout)

        alignment_info_var_1 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)

        if self.is_inseq:
            alignment_info_var_2 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_2.fa",
                self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)
        else:
            alignment_info_var_2 = copy.copy(alignment_info_var_1)


        if self.is_short_del_dup:
            with open(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_1.fa", 'w') as hout:
                print('>' + self.temp_key2 + '\n' + self.reference_segment_1, file = hout)

            with open(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_2.fa", 'w') as hout:
                print('>' + self.temp_key2 + '\n' + self.reference_segment_2, file = hout)

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
        if self.is_short_del_dup:
            supporting_reads = [rname for rname in supporting_reads if \
                (alignment_info_var_1[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres and \
                alignment_info_var_1[rname][0] >= alignment_info_ref_2[rname][0] + self.var_ref_margin_thres) or \
                (alignment_info_var_2[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres and \
                alignment_info_var_2[rname][0] >= alignment_info_ref_2[rname][0] + self.var_ref_margin_thres)]

        # filtering by mapping qualities
        supporting_reads = [rname for rname in supporting_reads if \
            self.readid2mapq[rname][0] is not None and self.readid2mapq[rname][0] >= self.var_read_min_mapq and \
            self.readid2mapq[rname][1] is not None and self.readid2mapq[rname][1] >= self.var_read_min_mapq]


        tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid = self.temp_key.split(',')
        print(f"{tchr1}\t{tpos1}\t{tdir1}\t{tchr2}\t{tpos2}\t{tdir2}\t{tinseq}\t{tid}\t{len(all_rnames)}\t{len(supporting_reads)}", file = self.hout_count)
        
        sread_info = { rname: alignment_info_var_1[rname] + alignment_info_var_2[rname] for rname in supporting_reads}
        for rname in supporting_reads:
            if self.is_short_del_dup:
                sinfo = '\t'.join([str(x) for x in alignment_info_var_1[rname] + alignment_info_var_2[rname] + alignment_info_ref_1[rname] + alignment_info_ref_2[rname]])
            else:
                sinfo = '\t'.join([str(x) for x in alignment_info_var_1[rname] + alignment_info_var_2[rname] + ["None"] * 12])
            print(f"{tchr1}\t{tpos1}\t{tdir1}\t{tchr2}\t{tpos2}\t{tdir2}\t{tinseq}\t{tid}\t{rname}\t{sinfo}", file = self.hout_ainfo)


    def count_alignment_and_print_sbnd(self):

        # close file hundle for writing long_read_seq.fa if it is open
        if self.temp_long_read_seq_file_h is not None: self.temp_long_read_seq_file_h.close()

        with open(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.variant_segment_1, file = hout)
        with open(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_1.fa", 'w') as hout:
            print('>' + self.temp_key2 + '\n' + self.reference_segment_1, file = hout)

        alignment_info_var_1 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".variant_seg_1.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)

        alignment_info_ref_1 = ssw_check(self.tmp_dir + '/' + self.temp_key2 + ".reference_seg_1.fa",
            self.tmp_dir + '/' + self.temp_key2 + ".long_read_seq.fa", self.use_ssw_lib)

        all_rnames = list(set(list(alignment_info_var_1.keys())))

        supporting_reads = [rname for rname in all_rnames if \
            (alignment_info_var_1[rname][0] > self.score_ratio_thres * len(self.variant_segment_1) and \
            alignment_info_var_1[rname][1] < self.start_pos_thres * len(self.variant_segment_1) and \
            alignment_info_var_1[rname][2] > self.end_pos_thres * len(self.variant_segment_1))]

        supporting_reads = [rname for rname in supporting_reads if \
            (alignment_info_var_1[rname][0] >= alignment_info_ref_1[rname][0] + self.var_ref_margin_thres)]

        # filtering by mapping qualities
        supporting_reads = [rname for rname in supporting_reads if \
            self.readid2mapq[rname][0] is not None and self.readid2mapq[rname][0] >= self.var_read_min_mapq]

        tchr, tpos, tdir, tcontig, tid = self.temp_key.split(',')
        print(f"{tchr}\t{tpos}\t{tdir}\t{tcontig}\t{tid}\t{len(all_rnames)}\t{len(supporting_reads)}", file = self.hout_count)

        sread_info = { rname: alignment_info_var_1[rname] for rname in supporting_reads}
        for rname in supporting_reads:
            sinfo = '\t'.join([str(x) for x in alignment_info_var_1[rname] ])
            print(f"{tchr}\t{tpos}\t{tdir}\t{tcontig}\t{tid}\t{rname}\t{sinfo}", file = self.hout_ainfo)



def count_sread_by_alignment(sv_file, alignment_file, output_count_file, output_alignment_info_file, reference_fasta, 
    sbnd_file = None, output_count_file_sbnd = None, output_alignment_info_file_sbnd = None,
    check_read_max_num = 500, var_read_min_mapq = 0, score_ratio_thres = 1.2, use_ssw_lib = False, sort_option = None, debug = False):

    gather_local_read_for_realignment(sv_file, alignment_file, output_count_file + ".tmp.local_read_for_realignment", reference_fasta,
        sbnd_file = sbnd_file, 
        output_file_sbnd = output_count_file_sbnd + ".tmp.local_read_for_realignment" if output_count_file_sbnd is not None else None,
        check_read_max_num = check_read_max_num,
        sort_option = sort_option)

    alignment_counter = Alignment_counter(output_count_file, output_alignment_info_file, reference_fasta,
        var_read_min_mapq, score_ratio_thres, use_ssw_lib, debug)

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

    if sbnd_file is None: return

    alignment_counter_sbnd = Alignment_counter(output_count_file_sbnd, output_alignment_info_file_sbnd, reference_fasta,
        var_read_min_mapq, score_ratio_thres, use_ssw_lib, debug)

    with open(output_count_file_sbnd + ".tmp.local_read_for_realignment", 'r') as hin:
        for line in hin:
            tkey, treadid, tmapq1, tmapq2, tseq = line.rstrip('\n').split('\t')
            tmapq2 = "None"
            if alignment_counter_sbnd.temp_key != tkey:
                if alignment_counter_sbnd.temp_key is not None:
                    alignment_counter_sbnd.count_alignment_and_print_sbnd()
                alignment_counter_sbnd.initialize(tkey, is_sbnd = True)

            alignment_counter_sbnd.add_long_read_seq(treadid, tseq, tmapq1, tmapq2)

        if alignment_counter_sbnd.temp_key is not None:
            alignment_counter_sbnd.count_alignment_and_print_sbnd()

    if not debug:
        os.remove(output_count_file_sbnd + ".tmp.local_read_for_realignment")

