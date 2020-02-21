#! /user/bin/env python

import sys
import pysam

def make_fasta_file(input_file, output_file, seq_id_file):

    
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout1, open(seq_id_file, 'w') as hout2:
        header = hin.readline()
        sid = 1
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if len(F[6]) < 100: continue
        
            key = ','.join(F[:6] + [str(len(F[6]))])

            print(">seq" + str(sid) + ',' + str(len(F[6])), file = hout1)
            print(F[6], file = hout1)
            sid = sid + 1
            print("seq" + str(sid) + '\t' + key, file = hout2)
            


def sam2bed_split(input_file, output_file):

    sam_file = pysam.AlignmentFile(input_file, 'r')
    hout = open(output_file, 'w')

    for read in sam_file.fetch():

        if read.is_unmapped: continue

        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary
        if read.is_secondary: continue

        query_name = read.query_name
        query_length = read.infer_read_length()
        query_strand = '-' if read.is_reverse else '+'

        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end

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

        query_pos_cur = query_start - 1 if query_strand == '+' else query_end
        query_pos_check = query_pos_cur
        query_exon_cur = query_pos_cur
        reference_pos_cur = reference_start - 1
        reference_exon_cur = reference_pos_cur
        reference_pos_check = reference_start - 1


        for cigar in cigartuples:
            if cigar[0] == 0:
                query_pos_cur = query_pos_cur + cigar[1] if query_strand == '+' else query_pos_cur - cigar[1]
                reference_pos_cur = reference_pos_cur + cigar[1] 
            elif cigar[0] == 1:
                query_pos_cur = query_pos_cur + cigar[1] if query_strand == '+' else query_pos_cur - cigar[1]
            elif cigar[0] == 2:
                reference_pos_cur = reference_pos_cur + cigar[1]
            elif cigar[0] == 3:
                exon_size = reference_pos_cur - reference_exon_cur 
                query_name2 = query_name + ',' + str(min(query_exon_cur, query_pos_cur) + 1) + ',' + \
                  str(max(query_exon_cur, query_pos_cur))
                print('\t'.join([reference_name, str(reference_exon_cur), str(reference_pos_cur), 
                                 query_name2, str(exon_size), query_strand]), file = hout)

                reference_pos_cur = reference_pos_cur + cigar[1]
                reference_exon_cur = reference_pos_cur
                query_exon_cur = query_pos_cur

        exon_size = reference_pos_cur - reference_exon_cur
        query_name2 = query_name + ',' + str(min(query_exon_cur, query_pos_cur) + 1) + ',' + \
        str(max(query_exon_cur, query_pos_cur))
        print('\t'.join([reference_name, str(reference_exon_cur), str(reference_pos_cur),
                         query_name2, str(exon_size), query_strand]), file = hout)

        if query_strand == '+' and query_end != query_pos_cur:
            print("query end inconsistent!!")
            sys.exit(1)
        if query_strand == '-' and query_start != query_pos_cur + 1:
            print("query end inconsistent!!")
            sys.exit(1)

    # read.qname = read.qname + ',' + str(query_start) + ',' + str(query_end) + ',' + query_strand
    hout.close()



def pp_proc_filt_exon(input_file, output_file):

    key2aln_info = {}
    key2insert_len = {}
    key2total_match = {}
    key2exon_match_num = {}

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            insert_key_info = F[3].split(',')
            insert_key = ','.join(insert_key_info[:7])
            insert_len = int(insert_key_info[6])
            aln_start, aln_end = int(insert_key_info[7]), int(insert_key_info[8])

            gene = F[9]
            exon_num = F[10]
            key = insert_key + '\t' + gene
            key2insert_len[key] = insert_len
            if key not in key2total_match: key2total_match[key] = 0
            if key not in key2exon_match_num: key2exon_match_num[key] = 0
            if key not in key2aln_info: key2aln_info[key] = []

            key2total_match[key] = key2total_match[key] + min(int(F[12]), aln_end - aln_start + 1)

            match_ratio1 = float(F[12]) / (int(F[8]) - int(F[7]))
            match_ratio2 = float(F[12]) / (int(F[2]) - int(F[1]))
            if min(match_ratio1, match_ratio2) >= 0.95: key2exon_match_num[key] = key2exon_match_num[key] + 1

            # if F[1] == "82247561":
            #     import pdb; pdb.set_trace()

            if int(F[7]) - int(F[1]) > 10:
                if F[5] == '+':
                    aln_start = min(aln_end, aln_start + (int(F[7]) - int(F[1])))
                else:
                    aln_end = max(aln_start, aln_end - (int(F[7]) - int(F[1])))
            
            if int(F[2]) - int(F[8]) > 10:
                if F[5] == '+':
                    aln_end = max(aln_start, aln_end - (int(F[2]) - int(F[8])))
                else:
                    aln_start = min(aln_end, aln_start + (int(F[2]) - int(F[8])))

            key2aln_info[key].append((aln_start, aln_end, exon_num, round(match_ratio1, 4)))

        for key in key2aln_info:

            total_match_ratio = float(key2total_match[key]) / key2insert_len[key]

            aln_bar = []
            for ainfo in sorted(key2aln_info[key]):
                aln_bar.append(','.join([str(ainfo[0]), str(ainfo[1]), str(ainfo[2]), str(ainfo[3])]))

            if total_match_ratio < 0.5: continue
            if key2exon_match_num[key] < 2: continue

            print(key + '\t' + str(round(total_match_ratio, 4)) + '\t' + str(key2exon_match_num[key]) + '\t' + ';'.join(aln_bar), file = hout)
     


def proc_rmsk_info(tkey, rmsk_info):

    _, total_len = tkey.split(',')
    total_len = int(total_len)
    alu_len, alu_count = 0, 0
    L1_len, L1_count, L1_dir = 0, 0, []
    SVA_len, SVA_count = 0, 0
    repeat_len, repeat_count = 0, 0
    polyAT_len = 0

    repeat_info = []
    for i in range(len(rmsk_info)):

        F = rmsk_info[i]

        if F[8] == 'C': F[8] = '-'
        if F[8] == '+':
            repeat_key = ','.join([F[5], F[6], F[8], F[9], F[10], F[11], F[12]])
        else:
            repeat_key = ','.join([F[5], F[6], F[8], F[9], F[10], F[12], F[13]])
        repeat_info.append(repeat_key)

        if F[10] == "LINE/L1":
            L1_len = L1_len + int(F[6]) - int(F[5]) + 1
            L1_count = L1_count + 1
            L1_dir.append(F[8])
        if F[10] == "SINE/Alu":
            alu_len = alu_len + int(F[6]) - int(F[5]) + 1
            alu_count = alu_count + 1
        if F[10] == "Retroposon/SVA":
            SVA_len = SVA_len + int(F[6]) - int(F[5]) + 1
            SVA_count = SVA_count + 1
        if F[9] == "(T)n" and i == 0:
            polyAT_len = polyAT_len + int(F[6])   
        if F[9] == "(A)n" and i == len(rmsk_info) - 1:
            polyAT_len = polyAT_len + total_len - int(F[5]) + 1


    repeat_class = "None"
    L1_ratio = min(1.0, float(L1_len) / (float(total_len) - float(polyAT_len)))
    Alu_ratio = min(1.0, float(alu_len) / (float(total_len) - float(polyAT_len)))
    SVA_ratio = min(1.0, float(SVA_len) / (float(total_len) - float(polyAT_len)))

    if L1_ratio >= 0.8:
        if L1_count == 1:
            repeat_class = "Simple_LINE1"
        elif L1_count == 2 and L1_dir[0] != L1_dir[1]:
            repeat_class = "Inverted_LINE1"
        else:
            repeat_class = "Other_LINE1"

    if Alu_ratio >= 0.8:
        repeat_class = "Alu"

    if SVA_ratio >= 0.8:
        repeat_class = "SVA"

    return([repeat_class, str(round(L1_ratio, 4)), str(round(Alu_ratio, 4)), str(round(SVA_ratio, 4)), repeat_info])


def summarize_rmsk(input_file, output_file):

    temp_key = ''
    temp_rmsk_info = []

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split()
            if len(F) <= 1 or F[0] in ["SW", "score"]: continue

            if temp_key != F[4]:
                if temp_key != '':

                    repeat_class, L1_ratio, Alu_ratio, SVA_ratio, repeat_info = proc_rmsk_info(temp_key, temp_rmsk_info)
                    print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + ',' + Alu_ratio + ',' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)

                temp_key = F[4]
                temp_rmsk_info = []
            

            temp_rmsk_info.append(F)
     
        if temp_key != '':
            repeat_class, L1_ratio, Alu_ratio, SVA_ratio, repeat_info = proc_rmsk_info(temp_key, temp_rmsk_info)
            print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + ',' + Alu_ratio + ',' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)



