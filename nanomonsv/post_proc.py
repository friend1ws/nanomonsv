#! /usr/bin/env python3

import sys, csv, itertools
import pysam, parasail

from nanomonsv.my_seq import reverse_complement
from nanomonsv.logger import get_logger

logger = get_logger(__name__)

class Sv(object):

    def __init__(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, sv_id, 
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl):

        self.chr1 = tchr1
        self.pos1 = tpos1
        self.dir1 = tdir1
        self.chr2 = tchr2
        self.pos2 = tpos2
        self.dir2 = tdir2
        self.inseq = tinseq
        self.sv_id = sv_id 
        self.total_read_tumor = total_read_tumor
        self.var_read_tumor = var_read_tumor
        self.total_read_ctrl = total_read_ctrl
        self.var_read_ctrl = var_read_ctrl
        self.filter = []


class Sv_filterer(object):

    def __init__(self, output_file, reference_fasta, is_control, simple_repeat_bed,  
        min_tumor_VAF = 0.05, min_indel_size = 50, bp_dist_margin = 30, validate_seg_len = 100):
        self.sv_list = []
        self.hout = open(output_file, 'w')
        self.bp_dist_margin = bp_dist_margin
        self.is_control = is_control
        self.reference_h = pysam.FastaFile(reference_fasta)
        self.min_tumor_VAF = min_tumor_VAF
        self.min_indel_size = min_indel_size
        self.bp_dist_margin = bp_dist_margin
        self.validate_seg_len = validate_seg_len
        self.simple_repeat_tb = pysam.TabixFile(simple_repeat_bed) if simple_repeat_bed is not None else None 
        self.header = None
        

    def __del__(self):
        self.hout.close()
        self.reference_h.close()


    def filter_close_both_breakpoints(self, sv1, sv2, filter_item = "Duplicate_with_close_SV"):

        # only apply when the first or the second sv is not insertion
        if len(sv1.inseq) >= self.min_indel_size or len(sv2.inseq) >= self.min_indel_size: return 

        if sv1.chr1 == sv2.chr1 and sv1.chr2 == sv2.chr2 and sv1.dir1 == sv2.dir1 and sv1.dir2 == sv2.dir2 and \
            abs(sv1.pos1 - sv2.pos1) < self.bp_dist_margin and abs(sv1.pos2 - sv2.pos2) < self.bp_dist_margin:

            if sv1.var_read_tumor is not None and sv2.var_read_tumor is not None:
                if sv1.var_read_tumor > sv2.var_read_tumor: 
                    sv1.filter.append(filter_item); return
                if sv2.var_read_tumor < sv1.var_read_tumor: 
                    sv2.filter.append(filter_item); return

            if len(sv1.inseq) < len(sv2.inseq):
                sv1.filter.append(filter_item); return
            if len(sv2.inseq) < len(sv1.inseq): 
                sv2.filter.append(filter_item); return

            if sv1.pos1 < sv2.pos1: 
                sv1.filter.append(filter_item); return
            if sv2.pos1 < sv1.pos1: 
                sv2.filter.append(filter_item); return

            if sv1.pos2 < sv2.pos2: 
                sv1.filter.append(filter_item); return
            if sv2.pos2 < sv1.pos2: 
                sv2.filter.append(filter_item) ; return

            sv2.filter.append(filter_item); return



    def filter_sv_insertion_match(self, sv, ins, filter_item = "Duplicate_with_insertion"):
    
        # only apply when the first sv is not insertion and the second sv is insertion type
        if len(sv.inseq) >= self.min_indel_size or len(ins.inseq) < self.min_indel_size: return 

        if not (sv.chr1 == ins.chr1 and abs(sv.pos1 - ins.pos1) <= self.bp_dist_margin) and \
            not (sv.chr2 == ins.chr2 and abs(sv.pos2 - ins.pos2) <= self.bp_dist_margin):
            return 

        ins_seg = self.reference_h.fetch(ins.chr1, max(ins.pos1 - self.validate_seg_len - self.bp_dist_margin - 1, 0), ins.pos1 - 1).upper()
        ins_seg = ins_seg + ins.inseq
        ins_seg = ins_seg + self.reference_h.fetch(ins.chr1, ins.pos2 - 1, ins.pos2 + self.validate_seg_len + self.bp_dist_margin - 1).upper()

        if sv.dir1 == '+':
            tseq = self.reference_h.fetch(sv.chr1, max(sv.pos1 - self.validate_seg_len - 1, 0), sv.pos1 - 1).upper()
        else:
            tseq = self.reference_h.fetch(sv.chr1, sv.pos1 - 1, sv.pos1 + self.validate_seg_len - 1).upper()
            tseq = reverse_complement(tseq)

        if sv.dir1 == '+':
            sv_seg = tseq + sv.inseq
        else:
            sv_seg = tseq + reverse_complement(sv.inseq)

        if sv.dir2 == '-':
            tseq = self.reference_h.fetch(sv.chr2, sv.pos2 - 1, sv.pos2 + self.validate_seg_len - 1).upper()
        else:
            tseq = self.reference_h.fetch(sv.chr2, max(sv.pos2 - self.validate_seg_len - 1, 0), sv.pos2 - 1).upper()
            tseq = reverse_complement(tseq)
 
        sv_seg = sv_seg + tseq

        user_matrix = parasail.matrix_create("ACGT", 2, -2)
        res = parasail.ssw(sv_seg, ins_seg, 3, 1, user_matrix)
        res_r = parasail.ssw(reverse_complement(sv_seg), ins_seg, 3, 1, user_matrix)

        if res.score1 > res_r.score1:
            match_ratio = float(res.score1) / (2 * (res.ref_end1 - res.ref_begin1 + 1))
            if match_ratio > 0.75:
                if res.read_begin1 < 0.1 * len(sv_seg) and res.read_end1 > 0.9 * len(sv_seg): 
                    sv.filter.append(filter_item)
                    return 
        else:
            match_ratio = float(res_r.score1) / (2 * (res_r.ref_end1 - res_r.ref_begin1 + 1))
            if match_ratio > 0.75:
                if res_r.read_begin1 < 0.1 * len(sv_seg) and res_r.read_end1 > 0.9 * len(sv_seg): 
                    sv.filter.append(filter_item)
                    return

        return


    def filter_dup_insertion(self, ins1, ins2, filter_item = "Duplicate_with_close_insertion"):

        # only apply when the first and the second sv is insertion type
        if len(ins1.inseq) < self.min_indel_size or len(ins2.inseq) < self.min_indel_size: return None

        bp_match = (ins1.chr1 == ins2.chr1 and abs(ins1.pos1 - ins2.pos1) <= 2 * self.bp_dist_margin) and \
            (ins1.chr2 == ins2.chr2 and abs(ins1.pos2 - ins2.pos2) <= 2 * self.bp_dist_margin)

        if bp_match == False: return None

        # insertion with shorter length will be filtered
        if len(ins1.inseq) < len(ins2.inseq):
            ins1.filter.append(filter_item)
            return
        else:
            ins2.filter.append(filter_item)


    def filter_sv_with_decoy(self, sv, filter_item = "SV_with_decoy"):

        if sv.chr1.endswith("decoy") or sv.chr2.endswith("decoy"): 
            sv.filter.append(filter_item)
            return


    def filter_small_sv(self, sv, filter_item = "Too_small_size"):

        if sv.chr1 == sv.chr2 and sv.dir1 == '+' and sv.dir2 == '-':
            sv_size = sv.pos2 - sv.pos1 + len(sv.inseq) - 1
            if sv_size < self.min_indel_size:
                sv.filter.append(filter_item)
                return

    def filter_low_vaf_sv(self, sv, filter_item = "Too_low_VAF"):

        if float(sv.var_read_tumor) / float(sv.total_read_tumor) < self.min_tumor_VAF:
            sv.filter.append(filter_item)

   
    def filter_simple_repeat(self, sv, filter_item = "Simple_repeat", simple_repeat_dist_margin = 30):

        if sv.chr1 == sv.chr2 and sv.dir1 == '+' and sv.dir2 == '-':
            sv_size = sv.pos2 - sv.pos1 + len(sv.inseq) - 1

            tabix_error_flag = False
            try:
                records = self.simple_repeat_tb.fetch(sv.chr1, max(sv.pos1 - simple_repeat_dist_margin + 1, 0), 
                    sv.pos1 + simple_repeat_dist_margin)
            except Exception as inst:
                print(f'{type(inst)}: {inst.args}', file = sys.stderr)
                tabix_error_flag = True

            if tabix_error_flag == False:
                for record_line in records:
                    record = record_line.split('\t')
                    if sv.pos1 >= int(record[1]) - simple_repeat_dist_margin and sv.pos2 <= int(record[2]) + simple_repeat_dist_margin:
                        sv.filter.append(filter_item)


            
    def add_sv(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, sv_id,
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl):
    
        sv = Sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, sv_id,
            total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl)
        self.sv_list.append(sv)


    def apply_filters(self):

        for sv in self.sv_list:
            self.filter_low_vaf_sv(sv)

        for sv in self.sv_list:
            self.filter_small_sv(sv)

        for sv in self.sv_list:
            self.filter_sv_with_decoy(sv)

        if self.simple_repeat_tb is not None:
            for sv in self.sv_list:
                self.filter_simple_repeat(sv)

        # logger.info("filter_close_both_breakpoints")
        for sv1, sv2 in itertools.combinations(self.sv_list, 2):
            if len(sv1.filter) > 0 or len(sv2.filter) > 0: continue
            self.filter_close_both_breakpoints(sv1, sv2)

        # logger.info("filter_sv_insertion_match")
        for sv1, sv2 in itertools.combinations(self.sv_list, 2):
            if len(sv1.filter) > 0 or len(sv2.filter) > 0: continue
            if len(sv1.inseq) < self.min_indel_size and len(sv2.inseq) >= self.min_indel_size: 
                self.filter_sv_insertion_match(sv1, sv2)
            elif len(sv2.inseq) < self.min_indel_size and len(sv1.inseq) >= self.min_indel_size:
                self.filter_sv_insertion_match(sv2, sv1)

        # logger.info("filter_dup_insertion")
        for ins1, ins2 in itertools.combinations(self.sv_list, 2):
            if len(ins1.filter) > 0 or len(ins2.filter) > 0: continue
            if len(ins1.inseq) >= self.min_indel_size and len(ins2.inseq) >= self.min_indel_size:
                self.filter_dup_insertion(ins1, ins2)


    def flush_sv_list(self):

        header = "Chr_1\tPos_1\tDir_1\tChr_2\tPos_2\tDir_2\tInserted_Seq\tSV_ID\tChecked_Read_Num_Tumor\tSupporting_Read_Num_Tumor" 
        if self.is_control: header = header + "\tChecked_Read_Num_Control\tSupporting_Read_Num_Control"
        header = header + '\t' + "Is_Filter"
        print(header, file = self.hout)
        for sv in self.sv_list:

            print_sv_line = f"{sv.chr1}\t{sv.pos1}\t{sv.dir1}\t{sv.chr2}\t{sv.pos2}\t{sv.dir2}\t{sv.inseq}\t{sv.sv_id}"
            print_sv_line = print_sv_line + f"\t{sv.total_read_tumor}\t{sv.var_read_tumor}"

            if self.is_control: print_sv_line = print_sv_line + f"\t{sv.total_read_ctrl}\t{sv.var_read_ctrl}"
            if len(sv.filter) == 0: 
                filter_print = "PASS"
            else:
                if "Duplicate_with_close_SV" in sv.filter: continue
                if "Duplicate_with_insertion" in sv.filter: continue
                if "Duplicate_with_close_insertion" in sv.filter: continue
                filter_print = ';'.join(list(set(sv.filter)))

            print(f"{print_sv_line}\t{filter_print}", file = self.hout)

    
def integrate_realignment_result(tumor_file, control_file, output_file, reference_fasta, simple_repeat_bed = None, min_indel_size = 50,
    min_tumor_variant_read_num = 3, min_tumor_VAF = 0.05, max_control_variant_read_num = 1, max_control_VAF = 0.03):

    is_control = True if control_file is not None else False
    svkey2control_info = {}
    if is_control:
        with open(control_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                svkey = (F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6], F[7])
                svkey2control_info[svkey] = (int(F[8]), int(F[9]))

    sv_filterer = Sv_filterer(output_file, reference_fasta, is_control, simple_repeat_bed, min_tumor_VAF, min_indel_size)
    with open(tumor_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6], F[7]
            tkey = (tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid)

            total_read_tumor, var_read_tumor = int(F[8]), int(F[9])
            total_read_ctrl, var_read_ctrl = None, None

            if tkey in svkey2control_info: total_read_ctrl, var_read_ctrl = svkey2control_info[tkey][0], svkey2control_info[tkey][1]

            # filtering by supporting read and variant frequencies
            if total_read_tumor == 0: continue
            if var_read_tumor < min_tumor_variant_read_num: continue

            if is_control:
                if total_read_ctrl is None or total_read_ctrl == 0: continue
                if var_read_ctrl > max_control_variant_read_num: continue
                if float(var_read_ctrl) / float(total_read_ctrl) > max_control_VAF: continue


            sv_filterer.add_sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, tid,
                total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl)

    sv_filterer.apply_filters()
    sv_filterer.flush_sv_list()
    del sv_filterer


def proc_sread_info_file(tumor_sread_info_file, sv_result_file, output_file, validate_sequence_length = 200):

    svkey = {}
    with open(sv_result_file, 'r') as hin:
        header = hin.readline()
        for line in hin:
            F = line.rstrip('\t').split('\t')
            svkey['\t'.join(F[:8])] = 1

    with open(tumor_sread_info_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tkey = '\t'.join(F[:8])
            if tkey not in svkey: continue

            score1, cstart1, cend1, sstart1, send1, strand1 = int(F[9]), int(F[10]), int(F[11]), int(F[12]), int(F[13]), F[14]
            score2, cstart2, cend2, sstart2, send2, strand2 = int(F[15]), int(F[16]), int(F[17]), int(F[18]), int(F[19]), F[20]

            tinseq = F[6]
            readid = F[8]
                
            if score1 >= score2:
                sstrand = strand1
                if sstrand == '+':
                    spos = (float(send1 - sstart1) / float(cend1 - cstart1)) * (validate_sequence_length - cstart1) + sstart1
                else:
                    spos = (float(sstart1 - send1) / float(cend1 - cstart1)) * (validate_sequence_length - cstart1) + send1
            else:
                sstrand = strand2
                if sstrand == '+':
                    spos = (float(send2 - sstart2) / float(cend2 - cstart2)) * (validate_sequence_length - cstart2) + sstart2 - len(tinseq) - 1
                else:
                    spos = (float(sstart2 - send2) / float(cend2 - cstart2)) * (validate_sequence_length - cstart2) + send2 + len(tinseq) + 1
            
            print(f"{tkey}\t{readid}\t{int(spos)}\t{sstrand}", file = hout)


def integrate_realignment_result_sbnd(tumor_sbnd_count_file, ctrl_sbnd_count_file, output_file,
    nonsbnd_result_file, refined_bp_sbnd_file,
    min_tumor_variant_read_num = 3, min_tumor_VAF = 0.05, max_control_variant_read_num = 1, max_control_VAF = 0.03):

    margin = 30
    nanomonsv_bp_list = {}
    with open(nonsbnd_result_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Chr_1": continue
            nanomonsv_bp_list[(F[0], int(F[1]) - margin, int(F[1]) + margin)] = 1
            nanomonsv_bp_list[(F[3], int(F[4]) - margin, int(F[4]) + margin)] = 1

    key2count_ctrl = {}
    if ctrl_sbnd_count_file is not None:
        with open(ctrl_sbnd_count_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                key = F[4]
                key2count_ctrl[key] = (F[5], F[6])

    key2contig = {}
    with open(refined_bp_sbnd_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[4]
            key2contig[key] = F[3]

    with open(tumor_sbnd_count_file, 'r') as hin, open(output_file, 'w') as hout:
        header = "Chr_1\tPos_1\tDir_1\tContig\tSV_ID\tChecked_Read_Num_Tumor\tSupporting_Read_Num_Tumor"
        if ctrl_sbnd_count_file is not None: header = header + "\tChecked_Read_Num_Control\tSupporting_Read_Num_Control"
        print(header, file = hout)

        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[4]
            """
            nanomonsv_flag = False
            for bp in nanomonsv_bp_list:
                if F[0] == bp[0] and int(F[1]) >= bp[1] and int(F[1]) <= bp[2]: nanomonsv_flag = True
            if nanomonsv_flag: continue
            """

            if ctrl_sbnd_count_file is not None:
                if key not in key2count_ctrl:
                    continue
                else:
                    ctrl_count = key2count_ctrl[key]

                if int(ctrl_count[0]) == 0: continue
                if int(ctrl_count[1]) > max_control_variant_read_num: continue
                if float(ctrl_count[1]) / float(ctrl_count[0]) > max_control_VAF: continue

            if int(F[6]) < min_tumor_variant_read_num: continue
            if float(F[6]) / float(F[5]) < min_tumor_VAF: continue

            if ctrl_sbnd_count_file is not None:
                print(f"{F[0]}\t{F[1]}\t{F[2]}\t{key2contig[key]}\t{F[4]}\t{F[5]}\t{F[6]}\t{ctrl_count[0]}\t{ctrl_count[1]}", file = hout)
            else:
                print(f"{F[0]}\t{F[1]}\t{F[2]}\t{key2contig[key]}\t{F[4]}\t{F[5]}\t{F[6]}", file = hout)
    

if __name__ == "__main__":

    import sys
    integrate_realignment_result(sys.argv[1], sys.argv[2], sys.argv[3])

