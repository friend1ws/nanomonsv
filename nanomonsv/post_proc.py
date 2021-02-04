#! /usr/bin/env python3


class Sv(object):

    def __init__(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, 
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl):

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


class Duplicate_remover(object):

    def __init__(self, output_file, is_control, bp_dist_margin = 20):
        self.sv_list = []
        self.hout = open(output_file, 'w')
        self.bp_dist_margin = bp_dist_margin
        self.is_control = is_control

    def __del__(self):
        self.hout.close()


    def check_duplicate(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq,
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl):

        for sv in self.sv_list:

            if tchr1 == sv.chr1 and tchr2 == sv.chr2 and tdir1 == sv.dir1 and tdir2 == sv.dir2 and \
                abs(tpos1 - sv.pos1) < self.bp_dist_margin and abs(tpos2 - sv.pos2) < self.bp_dist_margin:

                replace_flag = False
                if var_read_tumor > sv.var_read_tumor:
                    replace_flag = True
                elif var_read_tumor == sv.var_read_tumor:
                    if len(tinseq) < len(sv.inseq):
                        replace_flag = True
                    elif len(tinseq) == len(sv.inseq):
                        if tpos1 < sv.pos1 or tpos1 == sv.pos1 and tpos2 < sv.pos2:
                            replace_flag = True

                if replace_flag == True:
                    sv.pos1, sv.pos2, sv.inseq = tpos1, tpos2, tinseq
                    self.total_read_tumor, self.var_read_tumor = total_read_tumor, var_read_tumor
                    self.total_read_ctrl, self.var_read_ctrl = total_read_ctrl, var_read_ctrl

                return sv

        return None

    
    def add_sv(self, tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq,
        total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl):
                     
        sv = Sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, 
            total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl)
        self.sv_list.append(sv)


    def flush_sv_list(self):

        header = "Chr_1\tPos_1\tDir_1\tChr_2\tPos_2\tDir_2\tInserted_Seq\tChecked_Read_Num_Tumor\tSupporting_Read_Num_Tumor" 
        if self.is_control: header = header + "\tChecked_Read_Num_Control\tSupporting_Read_Num_Control"
        print(header, file = self.hout)
        for sv in self.sv_list:
            print_sv_line = f"{sv.chr1}\t{sv.pos1}\t{sv.dir1}\t{sv.chr2}\t{sv.pos2}\t{sv.dir2}\t{sv.inseq}"
            print_sv_line = print_sv_line + f"\t{sv.total_read_tumor}\t{sv.var_read_tumor}"
            if self.is_control: print_sv_line = print_sv_line + f"\t{sv.total_read_ctrl}\t{sv.var_read_ctrl}"
            print(print_sv_line, file = self.hout)

    
def integrate_realignment_result(tumor_file, control_file, output_file,
    min_tumor_variant_read_num = 3, min_tumor_VAF = 0.05, max_control_variant_read_num = 1, max_control_VAF = 0.03):

    is_control = True if control_file is not None else False
    svkey2control_info = {}
    if is_control:
        with open(control_file, 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                svkey = (F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6])
                svkey2control_info[svkey] = (int(F[7]), int(F[8]))

    duplicate_remover = Duplicate_remover(output_file, is_control)
    with open(tumor_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq = F[0], int(F[1]), F[2], F[3], int(F[4]), F[5], F[6]
            tkey = (tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq)

            total_read_tumor, var_read_tumor = int(F[7]), int(F[8])
            total_read_ctrl, var_read_ctrl = None, None

            if tkey in svkey2control_info: total_read_ctrl, var_read_ctrl = svkey2control_info[tkey][0], svkey2control_info[tkey][1]

            # filtering by supporting read and variant frequencies
            if total_read_tumor == 0: continue
            if var_read_tumor < min_tumor_variant_read_num: continue
            if float(var_read_tumor) / float(total_read_tumor) < min_tumor_VAF: continue

            if is_control:
                if total_read_ctrl is None or total_read_ctrl == 0: continue
                if var_read_ctrl > max_control_variant_read_num: continue
                if float(var_read_ctrl) / float(total_read_ctrl) > max_control_VAF: continue


            cret = duplicate_remover.check_duplicate(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq, 
                total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl)

            if cret is None:
                duplicate_remover.add_sv(tchr1, tpos1, tdir1, tchr2, tpos2, tdir2, tinseq,
                    total_read_tumor, var_read_tumor, total_read_ctrl, var_read_ctrl)

    duplicate_remover.flush_sv_list()
    del duplicate_remover


def proc_sread_info_file(tumor_sread_info_file, sv_result_file, output_file, validate_sequence_length = 200):

    svkey = {}
    with open(sv_result_file, 'r') as hin:
        header = hin.readline()
        for line in hin:
            F = line.rstrip('\t').split('\t')
            svkey['\t'.join(F[:7])] = 1

    with open(tumor_sread_info_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tkey = '\t'.join(F[:7])
            if tkey not in svkey: continue

            score1, cstart1, cend1, sstart1, send1, strand1 = int(F[8]), int(F[9]), int(F[10]), int(F[11]), int(F[12]), F[13]
            score2, cstart2, cend2, sstart2, send2, strand2 = int(F[14]), int(F[15]), int(F[16]), int(F[17]), int(F[18]), F[19]

            tinseq = F[6]
            readid = F[7]
                
            if score1 >= score2:
                sstrand = strand1
                if sstrand == '+':
                    spos = (float(send1 - sstart1) / float(cend1 - cstart1)) * (validate_sequence_length - cstart1) + sstart1
                else:
                    spos = (float(sstart1 - send1) / float(cend1 - cstart1)) * (validate_sequence_length - cstart1) + send1
            else:
                sstrand = strand2
                if sstrand == '+':
                    spos = (float(send1 - sstart1) / float(cend1 - cstart1)) * (validate_sequence_length - cstart1) + sstart1 - len(tinseq) - 1
                else:
                    spos = (float(sstart1 - send1) / float(cend1 - cstart1)) * (200 - cstart1) + send1 + len(tinseq) + 1
              
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
                key = f"{F[0]}\t{F[1]}\t{F[2]}"
                key2count_ctrl[key] = (F[4], F[5])

    key2contig = {}
    with open(refined_bp_sbnd_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = f"{F[0]}\t{F[1]}\t{F[2]}"
            key2contig[key] = F[3]

    with open(tumor_sbnd_count_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = f"{F[0]}\t{F[1]}\t{F[2]}"
            nanomonsv_flag = False
            for bp in nanomonsv_bp_list:
                if F[0] == bp[0] and int(F[1]) >= bp[1] and int(F[1]) <= bp[2]: nanomonsv_flag = True
                if F[3] == bp[0] and int(F[4]) >= bp[1] and int(F[4]) <= bp[2]: nanomonsv_flag = True

            if nanomonsv_flag: continue

            # import pdb; pdb.set_trace()
            if ctrl_sbnd_count_file is not None:
                if key not in key2count_ctrl:
                    continue
                else:
                    ctrl_count = key2count_ctrl[key]

                if int(ctrl_count[0]) == 0: continue
                if int(ctrl_count[1]) > max_control_variant_read_num: continue
                if float(ctrl_count[1]) / float(ctrl_count[0]) > max_control_VAF: continue

            if int(F[5]) < min_tumor_variant_read_num: continue
            if float(F[5]) / float(F[4]) < min_tumor_VAF: continue

            if ctrl_sbnd_count_file is not None:
                print(f"{key}\t{key2contig[key]}\t{F[4]}\t{F[5]}\t{ctrl_count[0]}\t{ctrl_count[1]}", file = hout)
            else:
                print(f"{key}\t{key2contig[key]}\t{F[4]}\t{F[5]}", file = hout)
    

if __name__ == "__main__":

    import sys
    integrate_realignment_result(sys.argv[1], sys.argv[2], sys.argv[3])

