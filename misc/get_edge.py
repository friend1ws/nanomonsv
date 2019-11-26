#! /usr/bin/env python

import sys, itertools

input_file = sys.argv[1]

consistent_pos_margin = 1.3
temp_rname = ""
temp_keys = []
temp_poses = []
temp_dirs = []
key2segment_size1 = {}
key2segment_size2 = {}

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        if F[1] != temp_rname:
            if temp_rname != "" and len(temp_keys) > 1:
                for i in range(0, len(temp_keys) - 1):

                    tchr1_1, tpos1_1, tdir1_1, tchr1_2, tpos1_2, tdir1_2, tinseq1 = temp_keys[i].split(',')
                    tchr2_1, tpos2_1, tdir2_1, tchr2_2, tpos2_2, tdir2_2, tinseq2 = temp_keys[i + 1].split(',')

                    jpos1, jpos2 = temp_poses[i], temp_poses[i + 1]
                    jdir1, jdir2 = temp_dirs[i], temp_dirs[i + 1]
            
                    if jdir1 == '+' and tinseq1 != "---": jpos1 = jpos1 + len(tinseq1)
                    if jdir2 == '+' and tinseq2 != "---": jpos2 = jpos2 - len(tinseq2)

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
                    # if spos1_2 == 25829678 or spos2_1 == 25829678:
                    #     import pdb; pdb.set_trace()

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
    
                    if consistent_flag == False:
                        if schr1_2 < schr2_1 or (schr1_2 == schr2_1 and spos1_2 > spos2_1):
                            schr1_1, schr2_2 = schr2_2, schr1_1
                            spos1_1, spos2_2 = spos2_2, spos1_1
                            sdir1_1, sdir2_2 = sdir2_2, sdir1_1
                            schr1_2, schr2_1 = schr2_1, schr1_2
                            spos1_2, spos2_1 = spos2_1, spos1_2
                            sdir1_2, sdir2_1 = sdir2_1, sdir1_2
                            sinseq1, sinseq2 = sinseq2, sinseq1

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

            temp_rname = F[1]
            temp_keys = []
            temp_poses = []
            temp_dirs = []

        temp_keys.append(F[0])
        temp_poses.append(int(F[2]))
        temp_dirs.append(F[3])


for key in key2segment_size1:
    segment_size1 = key2segment_size1[key]
    segment_size2 = ','.join(key2segment_size2[key])

    print('\t'.join([key, segment_size1, segment_size2]))

