#! /usr/bin/env python

import sys

def filter_support_reads(nanomonsv_results, support_read_file, output_file):
    pass_sv_ids = set()
    with open(nanomonsv_results, "r") as f:
        for i, line in enumerate(f):
            if i == 0: continue
            items = line.rstrip("\n").split("\t")
            sv_id = items[7]
            pass_sv_ids.add(sv_id)
    
    with open(support_read_file, "r") as f, open(output_file, "w") as w:
        for line in f:
            items = line.rstrip("\n").split("\t")
            sv_id = items[7]
            if sv_id in pass_sv_ids:
                print(line.rstrip("\n"), file=w)

def get_edge(sorted_support_read_file, output_file, consistent_pos_margin = 1.3):

    temp_rname = ""
    temp_keys = []
    temp_poses = []
    temp_dirs = []
    key2segment_size1 = {}
    key2segment_size2 = {}

    with open(sorted_support_read_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[-3] != temp_rname:
                if temp_rname != "" and len(temp_keys) > 1:
                    for i in range(0, len(temp_keys) - 1):
                        # Information of SV_1 and SV_2
                        tchr1_1, tpos1_1, tdir1_1, tchr1_2, tpos1_2, tdir1_2, tinseq1 = temp_keys[i].split(',')
                        tchr2_1, tpos2_1, tdir2_1, tchr2_2, tpos2_2, tdir2_2, tinseq2 = temp_keys[i + 1].split(',')

                        # Read position of SV breakpoint 1
                        jpos1, jpos2 = temp_poses[i], temp_poses[i + 1]
                        # Direction of realignment of reads against SV ref.
                        jdir1, jdir2 = temp_dirs[i], temp_dirs[i + 1]
                        
                        '''
                        Adjust read position of SV breakpoint
                        If realignment direction of SV_1 is "+", then read position is adjusted to read position of breakpoint 2.
                        If realignment direcion of SV_2 is "-", then read position is adjusted to read position of breakpoint 2.
                        '''
                        if jdir1 == '+' and tinseq1 != "---": jpos1 = jpos1 + len(tinseq1)
                        if jdir2 == '-' and tinseq2 != "---": jpos2 = jpos2 - len(tinseq2)
                        
                        # Modify SV order by realignment direction
                        if jdir1 == '+':
                            schr1_1, spos1_1, sdir1_1, schr1_2, spos1_2, sdir1_2 = \
                            tchr1_1, int(tpos1_1), tdir1_1, tchr1_2, int(tpos1_2), tdir1_2
                        else:
                            schr1_2, spos1_2, sdir1_2, schr1_1, spos1_1, sdir1_1 = \
                            tchr1_1, int(tpos1_1), tdir1_1, tchr1_2, int(tpos1_2), tdir1_2

                        if jdir2 == '+':
                            schr2_1, spos2_1, sdir2_1, schr2_2, spos2_2, sdir2_2 = \
                            tchr2_1, int(tpos2_1), tdir2_1, tchr2_2, int(tpos2_2), tdir2_2 
                        else:
                            schr2_2, spos2_2, sdir2_2, schr2_1, spos2_1, sdir2_1 = \
                            tchr2_1, int(tpos2_1), tdir2_1, tchr2_2, int(tpos2_2), tdir2_2

                        sinseq1, sinseq2 = tinseq1, tinseq2

                        # check positional consistency
                        segment_size = '---'
                        consistent_flag = False
                        segment_dir = '---'
                        if schr1_2 == schr2_1:
                            if sdir1_2 == '-' and sdir2_1 == '+': 
                                if spos2_1 - spos1_2 <= consistent_pos_margin * (jpos2 - jpos1) and \
                                spos2_1 - spos1_2 >= (1 / consistent_pos_margin) * (jpos2 - jpos1):
                                    consistent_flag = True
                                    segment_size = spos2_1 - spos1_2
                                    segment_dir = '+'
                            elif sdir1_2 == '+' and sdir2_1 == '-':
                                if spos1_2 - spos2_1 <= consistent_pos_margin * (jpos2 - jpos1) and \
                                spos1_2 - spos2_1 >= (1 / consistent_pos_margin) * (jpos2 - jpos1):
                                    consistent_flag = True
                                    segment_dir = '-'
                                    segment_size = spos1_2 - spos2_1

                        if consistent_flag == True and segment_dir == '-':
                            schr1_1, schr2_2 = schr2_2, schr1_1
                            spos1_1, spos2_2 = spos2_2, spos1_1
                            sdir1_1, sdir2_2 = sdir2_2, sdir1_1
                            schr1_2, schr2_1 = schr2_1, schr1_2
                            spos1_2, spos2_1 = spos2_1, spos1_2
                            sdir1_2, sdir2_1 = sdir2_1, sdir1_2
                            sinseq1, sinseq2 = sinseq2, sinseq1
        
                        if consistent_flag == False: continue
                        #    if schr1_2 < schr2_1 or (schr1_2 == schr2_1 and spos1_2 > spos2_1):
                        #        schr1_1, schr2_2 = schr2_2, schr1_1
                        #        spos1_1, spos2_2 = spos2_2, spos1_1
                        #        sdir1_1, sdir2_2 = sdir2_2, sdir1_1
                        #        schr1_2, schr2_1 = schr2_1, schr1_2
                        #        spos1_2, spos2_1 = spos2_1, spos1_2
                        #        sdir1_2, sdir2_1 = sdir2_1, sdir1_2
                        #        sinseq1, sinseq2 = sinseq2, sinseq1

                        key = '\t'.join([schr1_1, str(spos1_1), sdir1_1, sinseq1, schr1_2, str(spos1_2), sdir1_2,
                                        schr2_1, str(spos2_1), sdir2_1, sinseq2, schr2_2, str(spos2_2), sdir2_2, str(consistent_flag)])

                        if key not in key2segment_size1:
                            key2segment_size1[key], key2segment_size2[key] = str(segment_size), []

                        # key2segment_size1[key].append(str(segment_size))
                        key2segment_size2[key].append(str(jpos2 - jpos1))

                        """
                        print('\t'.join([schr1_1, str(spos1_1), sdir1_1, sinseq1, schr1_2, str(spos1_2), sdir1_2,
                                        schr2_1, str(spos2_1), sdir2_1, sinseq2, schr2_2, str(spos2_2), sdir2_2,
                                        str(consistent_flag), str(segment_size), str(jpos2 - jpos1)]))
                        """

                temp_rname = F[-3]
                temp_keys = []
                temp_poses = []
                temp_dirs = []

            temp_keys.append(",".join(F[0:7]))
            temp_poses.append(int(F[-2]))
            temp_dirs.append(F[-1])

    with open(output_file, "w") as w:
        for key in key2segment_size1:
            segment_size1 = key2segment_size1[key]
            segment_size2 = ','.join(key2segment_size2[key])

            print('\t'.join([key, segment_size1, segment_size2]), file=w)


def connect(edge_file, output_file):
    segments = []
    with open(edge_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            segments.append([F[:16] + ['+']])

    update = True
    while update:
        update = False
        for seg1 in range(len(segments) - 1):
            for seg2 in range(seg1 + 1, len(segments)):
                break1_1 = '\t'.join(segments[seg1][0][:3])
                break1_2 = '\t'.join(segments[seg1][0][4:7])
                break1_3 = '\t'.join(segments[seg1][-1][7:10])
                break1_4 = '\t'.join(segments[seg1][-1][11:14])

                break2_1 = '\t'.join(segments[seg2][0][:3])
                break2_2 = '\t'.join(segments[seg2][0][4:7])
                break2_3 = '\t'.join(segments[seg2][-1][7:10])
                break2_4 = '\t'.join(segments[seg2][-1][11:14])

                if break1_3 == break2_1 and break1_4 == break2_2:
                    sdir1, sdir2 = '+', '-'
                    update = True
                elif break1_3 == break2_4 and break1_4 == break2_3:
                    sdir1, sdir2 = '+', '+'
                    update = True
                elif break1_1 == break2_3 and break1_2 == break2_4:
                    sdir1, sdir2 = '-', '+'
                    update = True
                elif break1_1 == break2_2 and break1_2 == break2_1:
                    sdir1, sdir2 = '-', '-'
                    update = True

                if update: break

            if update: break


        if update:
            if sdir1 == sdir2:
                for j in range(len(segments[seg2])):
                    segments[seg2][j] = segments[seg2][j][11:14] + [segments[seg2][j][10]] + segments[seg2][j][7:10] + \
                                        segments[seg2][j][4:7] + [segments[seg2][j][3]] + segments[seg2][j][:3] + \
                                        segments[seg2][j][14:16] + ['+' if segments[seg2][j][16] == '-' else '-']

                    segments[seg2].reverse()

            if sdir1 == '+':
                segments[seg1] = segments[seg1] + segments[seg2]
            else:
                segments[seg1] = segments[seg2] + segments[seg1]

            del(segments[seg2])
            
    with open(output_file, "w") as w:
        for seg in segments:
            for i, s in enumerate(seg):
                if i == 0:
                    print("\t".join(s), end = "", file=w)
                else:
                    print("\t" + "\t".join(s[7:]), end = "", file=w)
                
                if i == len(seg) - 1:
                    print("", file=w)
    
def main():
    nanomonsv_results = sys.argv[1]
    support_read_file = sys.argv[2]
    filtered_support_read_file = sys.argv[3]
    output_edge_file = sys.argv[4]
    output_connect_file = sys.argv[5]
    
    filter_support_reads(nanomonsv_results, support_read_file, filtered_support_read_file)
    get_edge(filtered_support_read_file, output_edge_file)
    connect(output_edge_file, output_connect_file)

if __name__ == "__main__":
    main()  


"""
8       130933956       +       ---     8       130871365       -       8       130872760       +       ---     8       131401249       +       True    1395 
"""



