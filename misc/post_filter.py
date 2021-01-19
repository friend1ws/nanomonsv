#! /usr/bin/env python3

import csv, itertools
import pysam, parasail

from nanomonsv.my_seq import reverse_complement
from nanomonsv.logger import get_logger

logger = get_logger(__name__)

class Sv(object):

    def __init__(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, 
        total_read_tumor = None, var_read_tumor = None, total_read_ctrl = None, var_read_ctrl = None, print_line = None,
        insert_type = None):

        self.chr1 = tchr1
        self.pos1 = tpos1
        self.dir1 = tdir1
        self.chr2 = tchr2
        self.pos2 = tpos2
        self.dir2 = tdir2
        self.inseq = tinseq
        self.total_read_tumor = total_read_tumor
        self.var_read_tumor = var_read_tumor
        self.total_read_ctrl = total_read_ctrl
        self.var_read_ctrl = var_read_ctrl
        self.insert_type = insert_type
        self.print_line = print_line
        self.filter = []


class Duplicate_remover(object):

    def __init__(self, output_file, reference_fasta, simple_repeat_bed = None, 
        bp_dist_margin = 30, validate_seg_len = 100, simple_repeat_dist_margin = 50):
        self.sv_list = []
        self.hout = open(output_file, 'w')
        self.reference_h = pysam.FastaFile(reference_fasta)
        self.simple_repeat_tb = pysam.TabixFile(simple_repeat_bed) if simple_repeat_bed is not None else None
        self.bp_dist_margin = bp_dist_margin
        self.validate_seg_len = validate_seg_len
        self.simple_repeat_dist_margin = simple_repeat_dist_margin
        self.header = None

    def __del__(self):
        self.hout.close()
        self.reference_h.close()
        if self.simple_repeat_tb is not None: self.simple_repeat_tb.close()

    def add_sv(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq,
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl, print_line):
                     
        sv = Sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, 
            total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl, print_line)
        self.sv_list.append(sv)


    def filter_close_both_breakpoints(self, sv1, sv2, filter_item = "Duplicate_with_close_SV"):

        if sv1.chr1 == sv2.chr1 and sv1.chr2 == sv2.chr2 and sv1.dir1 == sv2.dir1 and sv1.dir2 == sv2.dir2 and \
            abs(sv1.pos1 - sv2.pos1) < self.bp_dist_margin and abs(sv1.pos2 - sv2.pos2) < self.bp_dist_margin:

            second_flag = False
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
        if len(sv.inseq) >= 100 or len(ins.inseq) < 100: return 

        if not (sv.chr1 == ins.chr1 and abs(sv.pos1 - ins.pos1) <= self.bp_dist_margin) and \
            not (sv.chr2 == ins.chr2 and abs(sv.pos2 - ins.pos2) <= self.bp_dist_margin):
            return 

        ins_seg = self.reference_h.fetch(ins.chr1, max(ins.pos1 - self.validate_seg_len - self.bp_dist_margin - 1, 0), ins.pos1 - 1)
        ins_seg = ins_seg + ins.inseq
        ins_seg = ins_seg + self.reference_h.fetch(ins.chr1, ins.pos2 - 1, ins.pos2 + self.validate_seg_len + self.bp_dist_margin - 1)

        if sv.dir1 == '+':
            tseq = self.reference_h.fetch(sv.chr1, max(sv.pos1 - self.validate_seg_len - 1, 0), sv.pos1 - 1)
        else:
            tseq = self.reference_h.fetch(sv.chr1, sv.pos1 - 1, sv.pos1 + self.validate_seg_len - 1) 
            tseq = reverse_complement(tseq)

        if sv.dir1 == '+':
            sv_seg = tseq + sv.inseq
        else:
            sv_seg = tseq + reverse_complement(sv.inseq)

        if sv.dir2 == '-':
            tseq = self.reference_h.fetch(sv.chr2, sv.pos2 - 1, sv.pos2 + self.validate_seg_len - 1)
        else:
            tseq = self.reference_h.fetch(sv.chr2, max(sv.pos2 - self.validate_seg_len - 1, 0), sv.pos2 - 1)
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

        # only apply when the first sv is not insertion and the second sv is insertion type
        if len(ins1.inseq) < 100 or len(ins2.inseq) < 100: return None

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


    def filter_small_sv(self, sv, size_thres = 100, filter_item = "Too_small_SV"):

        if sv.chr1 == sv.chr2 and sv.dir1 == '+' and sv.dir2 == '-':
            sv_size = sv.pos2 - sv.pos1 + len(sv.inseq) - 1
            if sv_size < size_thres:
                sv.filter.append(filter_item)
                return

    def filter_indel_in_simple_repeat(self, sv, rescure_annotated_insertions = True, filter_item = "Simple_repeat"):

        if sv.chr1 == sv.chr2 and sv.dir1 == '+' and sv.dir2 == '-':
            sv_size = sv.pos2 - sv.pos1 + len(sv.inseq) - 1
 
            if sv.insert_type is not None and rescure_annotated_insertions: return

            tabix_error_flag = False
            try:
                records = self.simple_repeat_tb.fetch(sv.chr1, max(sv.pos1 - self.simple_repeat_dist_margin + 1, 0), 
                    sv.pos1 + self.simple_repeat_dist_margin)
            except Exception as inst:
                logger.warning("%s: %s" % (type(inst), inst.args))
                tabix_error_flag = True

            if tabix_error_flag == False:
                for record_line in records:
                    record = record_line.split('\t')
                    if sv.pos1 >= int(record[1]) - self.simple_repeat_dist_margin and \
                        int(sv.pos2) <= int(record[2]) + self.simple_repeat_dist_margin:
                        sv.filter.append(filter_item); return


    def apply_filters(self):

        # logger.info("filter_close_both_breakpoints")
        for sv1, sv2 in itertools.combinations(self.sv_list, 2):
            self.filter_close_both_breakpoints(sv1, sv2)

        # logger.info("filter_sv_insertion_match")
        for sv1, sv2 in itertools.combinations(self.sv_list, 2):

            if len(sv1.inseq) < 100 and len(sv2.inseq) >= 100: 
                self.filter_sv_insertion_match(sv1, sv2)
            elif len(sv2.inseq) < 100 and len(sv1.inseq) >= 100:
                self.filter_sv_insertion_match(sv2, sv1)

        # logger.info("filter_dup_insertion")
        for ins1, ins2 in itertools.combinations(self.sv_list, 2):

            if len(ins1.inseq) >= 100 and len(ins2.inseq) >= 100:
                self.filter_dup_insertion(ins1, ins2)

        for sv in self.sv_list:
            self.filter_sv_with_decoy(sv)

        for sv in self.sv_list:
            self.filter_small_sv(sv)

        if self.simple_repeat_tb is not None:
            for sv in self.sv_list:
                self.filter_indel_in_simple_repeat(sv)
    

    def flush_sv_list(self):

        print(self.header + "\tIs_filter", file = self.hout)
        for sv in self.sv_list:
            if len(sv.filter) == 0: 
                filter_print = "PASS"
            else:
                filter_print = ';'.join(sv.filter)

            print(f"{sv.print_line}\t{filter_print}", file = self.hout)


def post_filter_main(args):

    duplicate_remover = Duplicate_remover(args.output_file, args.reference_fasta, simple_repeat_bed = args.simple_repeat_bed)

    with open(args.sv_list_file, 'r') as hin:
        dreader = csv.DictReader(hin, delimiter = '\t')
        header = dreader.fieldnames
        duplicate_remover.header = '\t'.join(header)
        if "Checked_Read_Num_Control" in header and "Supporting_Read_Num_Control" in header: is_control = True
        if "Insert_Type" in header: is_insert_type = True

        for F in dreader:
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F["Chr_1"], int(F["Pos_1"]), F["Dir_1"], F["Chr_2"], int(F["Pos_2"]), F["Dir_2"], F["Inserted_Seq"]

            total_read_tumor, var_read_tumor = int(F["Checked_Read_Num_Tumor"]), int(F["Supporting_Read_Num_Tumor"])
            if is_control:
                total_read_ctrl, var_read_ctrl = int(F["Checked_Read_Num_Control"]), int(F["Supporting_Read_Num_Control"])
            else:
                total_read_ctrl, var_read_ctrl = None, None

            """
            # filtering by supporting read and variant frequencies
            if total_read_tumor == 0: continue
            if var_read_tumor < min_tumor_variant_read_num: continue
            if float(var_read_tumor) / float(total_read_tumor) < min_tumor_VAF: continue

            if is_control:
                if total_read_ctrl is None or total_read_ctrl == 0: continue
                if var_read_ctrl > max_control_variant_read_num: continue
                if float(var_read_ctrl) / float(total_read_ctrl) > max_control_VAF: continue
            """

            duplicate_remover.add_sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq,
                total_read_tumor = total_read_tumor, var_read_tumor = var_read_tumor, 
                total_read_ctrl = total_read_ctrl, var_read_ctrl = var_read_ctrl, 
                print_line = '\t'.join(F.values()))


    duplicate_remover.apply_filters()
    duplicate_remover.flush_sv_list()
    del duplicate_remover


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(prog = "nanomonsv_post_filter",
        description = "Post filtering script for the result of nanomonsv")

    parser.add_argument("sv_list_file", type = str,
                        help = "Path to the nanomonsv result file")

    parser.add_argument("output_file", type = str,
                        help = "Path to the output file")

    parser.add_argument("reference_fasta", metavar = "reference.fa", type = str,
                        help = "Path to the reference genome sequence")

    parser.add_argument("--simple_repeat_bed", metavar = "simpleRepeat.bed.gz", type = str, default = None,
                        help = "Path to the tabix indexed simple repeat bed file")

    args = parser.parse_args()

    post_filter_main(args)


