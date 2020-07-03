#! /user/bin/env python

import sys, re, pkg_resources
import pysam

from .swalign import *
from .my_seq import get_seq, reverse_complement

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
            print("seq" + str(sid) + ',' + str(len(F[6])) + '\t' + key, file = hout2)
            sid = sid + 1
 


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



def pp_proc_filt_exon(input_file, seq_list, output_file):

    key2aln_info = {}
    key2insert_len = {}
    key2total_match = {}
    key2exon_match_num = {}

    sid2skey = {}
    with open(seq_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            sid2skey[F[0]] = F[1]

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            insert_key_info = F[3].split(',')
            insert_key = ','.join(insert_key_info[:2])
            insert_len = int(insert_key_info[1])
            aln_start, aln_end = int(insert_key_info[2]), int(insert_key_info[3])

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
                    print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + '\t' + Alu_ratio + '\t' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)

                temp_key = F[4]
                temp_rmsk_info = []
            

            temp_rmsk_info.append(F)
     
        if temp_key != '':
            repeat_class, L1_ratio, Alu_ratio, SVA_ratio, repeat_info = proc_rmsk_info(temp_key, temp_rmsk_info)
            print(temp_key + '\t' + repeat_class + '\t' + L1_ratio + '\t' + Alu_ratio + '\t' + SVA_ratio + '\t' + ';'.join(repeat_info), file = hout)



def check_tsd_polyAT(input_file, seq_list, reference, output_file):

    match = 1
    mismatch = -8
    scoring = NucleotideScoringMatrix(match, mismatch)
    sw = LocalAlignment(scoring, gap_penalty=-8, gap_extension_penalty=-8)

    sid2skey = {}
    with open(seq_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            sid2skey[F[0]] = F[1]

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            sid = line.rstrip('\n').lstrip('>')
            seq = hin.readline().rstrip('\n')
            skey = sid2skey[sid]

            tchr, tpos1, _, _, tpos2, _, tlen = skey.split(',')

            match1 = re.search(r'^([ACGT]{1,20}?)TTTTTTTTTT', seq)
            match2 = re.search(r'AAAAAAAAAA([ACGT]{1,20}?)$', seq)
            tsd1 = None
            tsd2 = None
            is_polyA = False
            is_polyT = False

            if seq.startswith('TTTTTTTTTT') or match1 is not None: is_polyT = True
            if seq.endswith('AAAAAAAAAA') or match2 is not None: is_polyA = True


            tsd_cand1 = seq[:20]
            local_seq1 = get_seq(reference, tchr, int(tpos2), int(tpos2) + len(tsd_cand1))
            alignment1 = sw.align(tsd_cand1, local_seq1)
            if alignment1.q_pos <= 2 and alignment1.matches >= 5 and alignment1.identity >= 0.8: 
                tsd1 = alignment1.query[(alignment1.q_pos):(alignment1.q_pos + alignment1.matches)]  


            tsd_cand2 = reverse_complement(seq[-20:])
            local_seq2 = reverse_complement(get_seq(reference, tchr, int(tpos1) - len(tsd_cand2), int(tpos1)))
            alignment2 = sw.align(tsd_cand2, local_seq2)
            if alignment2.q_pos <= 2 and alignment2.matches >= 5 and alignment2.identity >= 0.8: 
                tsd2 = reverse_complement(alignment2.query[(alignment2.q_pos):(alignment2.q_pos + alignment2.matches)])

            if tsd1 is not None and tsd2 is None: 
                tsd = tsd1
            elif tsd1 is None and tsd2 is not None:
                tsd = tsd2
            elif tsd1 is not None and tsd2 is not None:
                tsd = tsd1 if len(tsd1) >= len(tsd2) else tsd2
            else:
                tsd = None

            polyAT = None
            if is_polyT: polyAT = "polyT"
            if is_polyA: polyAT = "polyA"

            print(sid + '\t' + str(polyAT) + '\t' + str(tsd) + '\t' + tsd_cand1 + '\t' + local_seq1 + '\t' + tsd_cand2 + '\t' + local_seq2, file = hout) 



def summarize_bwa_alignment(input_sam, seq_list, output_file):

    samfile = pysam.AlignmentFile(input_sam, 'r')

    inserted_positions = []
    with open(seq_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr, tstart, _, _, tend, _, _ = F[1].split(',')
            inserted_positions.append((tchr, int(tstart), int(tend)))


    samfile = pysam.AlignmentFile(input_sam, 'r')
    hout = open(output_file, 'w') 
    for read in samfile.fetch():

        if read.is_secondary or read.is_supplementary: continue
        mapping_quality = read.mapping_quality
        if mapping_quality < 30: continue
 
        query_name = read.query_name
        key = query_name

        query_strand = '-' if read.is_reverse else '+'
        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end
        mapping_quality = read.mapping_quality
        query_length = read.infer_read_length()

        if query_strand == '+':
            query_start = read.query_alignment_start + 1 
            query_end = read.query_alignment_end
        else:
            query_start = query_length - read.query_alignment_end + 1
            query_end = query_length - read.query_alignment_start

        inserted_pos = "---"
        for inserted_position in inserted_positions:
            ichr, istart, iend = inserted_position
            if reference_name == ichr and reference_end > istart - 5000 and reference_start < iend + 5000:
                inserted_pos = ichr + ',' + str(istart) + ',' + str(iend)

        print(query_name + '\t' + ','.join([str(query_start), str(query_end), query_strand, str(mapping_quality), reference_name, str(reference_start), str(reference_end)]) + '\t' + inserted_pos, file = hout)

    samfile.close()
    hout.close()



def summarize_bwa_alignment2(input_sam, seq_list, output_file):

    samfile = pysam.AlignmentFile(input_sam, 'r')


    inserted_positions = []
    with open(seq_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr, tstart, _, _, tend, _, _ = F[1].split(',')
            inserted_positions.append((tchr, int(tstart), int(tend)))


    samfile = pysam.AlignmentFile(input_sam, 'r')

    query2align_supp = {}
    query2align_primary = {}
    query2align_ins = {}

    hout = open(output_file, 'w') 
    for read in samfile.fetch():

        if read.is_unmapped: continue
        if read.is_secondary: continue

        query_name = read.query_name
        key = query_name

        query_strand = '-' if read.is_reverse else '+'
        reference_name = read.reference_name
        reference_start = read.reference_start + 1
        reference_end = read.reference_end
        mapping_quality = read.mapping_quality
        query_length = read.infer_read_length()

        cigartuples = read.cigartuples
        left_hard_clipping_size, right_hard_clipping_size = 0, 0
        if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
        if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]

        if not read.is_supplementary:
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

        align_info = ','.join([str(query_start), str(query_end), query_strand, str(mapping_quality), reference_name, str(reference_start), str(reference_end)])


        if read.is_supplementary:
            if query_name not in query2align_supp: query2align_supp[query_name] = []
            query2align_supp[query_name].append(align_info + ",s")
        else:
            if read.mapping_quality >= 30:
                query2align_primary[query_name] = align_info + ",p"

                inserted_pos = "---"
                for inserted_position in inserted_positions:
                    ichr, istart, iend = inserted_position
                    if reference_name == ichr and reference_end > istart - 5000 and reference_start < iend + 5000:
                        inserted_pos = ichr + ',' + str(istart) + ',' + str(iend)
                query2align_ins[query_name] = inserted_pos

    samfile.close()

    hout = open(output_file, 'w')
    for query in sorted(query2align_primary):
        align_primary = query2align_primary[query]
        _, _, _, _, pchr, pstart, pend, _ = align_primary.split(',')
        align_ins = query2align_ins[query]
        align_info_final = align_primary
        if query in query2align_supp:
            for elm in query2align_supp[query]:
                _, _, _, _, schr, sstart, send, _ = elm.split(',')
                if int(pstart) - 5000 <= int(send) and int(pend) + 5000 >= int(sstart):
                    align_info_final = align_info_final + ';' + elm
        print("%s\t%s\t%s" % (query, align_info_final, align_ins), file = hout)


    hout.close()




def organize_info(rmsk_file, alignment_file, tsd_file, seq_list, output_file, genome_id):

    key2rmsk = {}
    with open(rmsk_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key2rmsk[F[0]] = [F[1], F[2], F[3], F[4], F[5]]

    key2alignment = {}
    with open(alignment_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key2alignment[F[0]] = [F[1], F[2]]


    key2tsd_polyAT = {}
    with open(tsd_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key2tsd_polyAT[F[0]] = [F[1], F[2]]

    keys = list(set(list(key2rmsk) + list(key2alignment) + list(key2tsd_polyAT)))

    if genome_id == "hg38":
        line1_tb = pysam.TabixFile(pkg_resources.resource_filename("nanomonsv", "data/LINE1.hg38.bed.gz"))
    else:
        line1_tb = pysam.TabixFile(pkg_resources.resource_filename("nanomonsv", "data/LINE1.hg19.bed.gz"))

    with open(output_file, 'w') as hout:

        for key in keys:

            repeat_type, line1_ratio, alu_ratio, sva_ratio, rmsk_info = key2rmsk[key] if key in key2rmsk else ["None", 0.0, 0.0, 0.0, "---"]
            alignment_infos, inserted_pos = key2alignment[key] if key in key2alignment else ["---", "---"]
            is_polyAT, tsd = key2tsd_polyAT[key] if key in key2tsd_polyAT else ["---", "---"]
            
            line1_ratio = float(line1_ratio)
            line1_info = None
            overlap_ratio = None
            transduction_class = None

            if repeat_type.endswith("LINE1"): transduction_class = "Solo"

            if alignment_infos == "---" or is_polyAT == None: 
                print('\t'.join([key, repeat_type, str(line1_ratio), str(alu_ratio), str(sva_ratio), 
                                 rmsk_info, alignment_infos, inserted_pos, is_polyAT, tsd,
                                 str(line1_info), str(transduction_class)]), file = hout)
                continue


            achr, astart, aend, astrand = None, None, None, None
            aqstart, aqend = None, None
            for alignment_info in alignment_infos.split(';'):
                tqstart, tqend, tastrand, _, tachr, tastart, taend, _ = alignment_info.split(',')
                if achr is None:
                    achr, astart, aend, astrand, aqstart, aqend = tachr, int(tastart), int(taend), tastrand, int(tqstart), int(tqend)
                else:
                    if (is_polyAT == "polyT" and int(tqstart) < aqstart) or (is_polyAT == "polyA" and int(tqend) > aqend):
                        achr, astart, aend, astrand, aqstart, aqend = tachr, int(tastart), int(taend), tastrand, int(tqstart), int(tqend)

            if not achr.startswith("chr"): achr = "chr" + achr

            source_dir = '+' if (astrand == '+' and is_polyAT == "polyT") or (astrand == '-' and is_polyAT == "polyA") else '-'

            tabix_error_flag = 0
            if source_dir == '+':
                cur_L1_pos = float("Inf")
                is_rmsk, is_1000g, is_gnomad = False, False, False
                try:
                    records = line1_tb.fetch(achr, int(astart) - 50, int(aend) + 5000)
                except Exception as inst:
                    tabix_error_flag = 1
                if tabix_error_flag == 0:
                    for record in records:
                        FF = record.split('\t')
                        if FF[5] == '+': continue
                        if int(FF[1]) > cur_L1_pos: continue
                        line1_info = FF[3]

                        if "umary_LINE1" in line1_info:
                            if is_rmsk: continue
                            is_1000g = True
                        elif "gnomAD-SV" in line1_info:
                            if is_1000g or is_rmsk: continue
                            is_gnomad = True
                        else:
                            is_rmsk = True

                        if line1_ratio < 0.01:
                            transduction_class = "Orphan"
                        else:
                            transduction_class = "Partnered"
                        cur_L1_pos = int(FF[1])

            elif source_dir == '-':
                cur_L1_pos = -float("Inf")
                is_rmsk, is_1000g, is_gnomad = False, False, False
                try:
                    records = line1_tb.fetch(achr, int(astart) - 5000, int(aend) + 50)
                except Exception as inst:
                    tabix_error_flag = 1
                if tabix_error_flag == 0:
                    for record in records:
                        FF = record.split('\t')
                        if FF[5] == '-': continue
                        if int(FF[1]) < cur_L1_pos: continue 
                        line1_info = FF[3]

                        if "umary_LINE1" in line1_info:                  
                            if is_rmsk: continue
                            is_1000g = True
                        elif "gnomAD-SV" in line1_info:
                            if is_1000g or is_rmsk: continue
                            is_gnomad = True
                        else:
                            is_rmsk = True

                        if line1_ratio < 0.01:
                            transduction_class = "Orphan"
                        else:
                            transduction_class = "Partnered"
                        cur_L1_pos = int(FF[1])


            print('\t'.join([key, repeat_type, str(line1_ratio), str(alu_ratio), str(sva_ratio),
                             rmsk_info, alignment_infos, inserted_pos, is_polyAT, tsd, 
                             str(line1_info), str(transduction_class)]), file = hout)
 


def annotate_sv_file(sv_file, source_file, ppseudo_file, seq_list, output_file):

    sid2skey = {}
    with open(seq_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            sid2skey[F[0]] = F[1]

    skey2source_info = {}
    with open(source_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            skey = sid2skey[F[0]]
            skey2source_info[skey] = F[1:]

    skey2ppseudo_info = {}
    with open(ppseudo_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            skey = sid2skey[F[0]]
            if skey in skey2ppseudo_info:
                tmatch_ratio = float(skey2ppseudo_info[skey][1])
                if float(F[2]) > tmatch_ratio:
                    skey2ppseudo_info[skey] = F[1:4]
            else:
                skey2ppseudo_info[skey] = F[1:4]


    with open(sv_file, 'r') as hin, open(output_file, 'w') as hout:
        header = hin.readline().rstrip('\n')
        print(header + '\t' + '\t'.join(["Insert_Type", "Is_Inversion", "L1_Raito", "Alu_Ratio", "SV_Ratio", "RMSK_Info",
                                         "Alignment_Info", "Inserted_Pos", "Is_PolyA_T", "Target_Site_Duplication", "L1_Source_Info",
                                         "PSD_Gene", "PSD_Overlap_Ratio", "PDS_Exon_Num"]), file = hout) 
        for line in hin:
            F = line.rstrip('\n').split('\t')
            skey = ','.join(F[:6] + [str(len(F[6]))])

            if skey not in skey2source_info:
                print('\t'.join(F) + '\t' + "None" + '\t' + '\t'.join(["---"] * 13), file = hout)
                continue
        
            source_info = skey2source_info[skey]
            ppseudo_line = '\t'.join(skey2ppseudo_info[skey]) if skey in skey2ppseudo_info else "---\t---\t---"

            insert_type = "---"
            if source_info[10] == "Orphan":
                insert_type = "Orphan_L1"
            elif source_info[10] == "Partnered":
                insert_type = "Partnered_L1"
            elif source_info[10] == "Solo":
                insert_type = "Solo_L1"
            elif source_info[0] == "Alu":
                insert_type = "Alu"
            elif source_info[0] == "SVA":
                insert_type = "SVA"
            elif skey in skey2ppseudo_info:
                insert_type = "PSD"

            is_inversion = "NA"
            if source_info[0].startswith("Simple"):
                is_inversion = "Simple"
            elif source_info[0].startswith("Inverted"):
                is_inversion = "Inverted"
            elif source_info[0].startswith("Other"):
                is_inversion = "Other"

            print('\t'.join(F) + '\t' + insert_type + '\t' + is_inversion + '\t' + \
                  '\t'.join(source_info[1:10]) + '\t' + ppseudo_line, file = hout)



