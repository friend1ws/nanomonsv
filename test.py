#! /usr/bin/env python

import sys
import pysam

bamfile = pysam.AlignmentFile(sys.argv[1], 'r')

temp_key = ''
bp1_info = None
bp2_info = None
for read in bamfile.fetch():

    if temp_key != read.query_name:

        if temp_key != '':
            print(temp_key, bp1_info, bp2_info)

        temp_key = read.query_name
        # if temp_key == "15,25831071,25831113,+,15,27447753,27447796,-":
        #     import pdb; pdb.set_trace()

        bp1_info = None
        bp2_info = None
        bchr1, bstart1, bend1, bdir1, bchr2, bstart2, bend2, bdir2 = temp_key.split(',')
        bstart1, bend1, bstart2, bend2 = int(bstart1), int(bend1), int(bstart2), int(bend2)
        
    
    query_strand = '-' if read.is_reverse else '+'
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


    reference_name = read.reference_name
    reference_start = read.reference_start + 1
    reference_end = read.reference_end

    # check bp1
    if reference_name == bchr1:
        if bdir1 == '+' and reference_end >= bstart1 and reference_end <= bend1:
            bp1_info = (reference_name, reference_end, '+', query_start, query_end)
        if bdir1 == '-' and reference_start >= bstart1 and reference_start <= bend1:
            bp1_info = (reference_name, reference_start, '-', query_start, query_end)


    if reference_name == bchr2:
        if bdir2 == '+' and reference_end >= bstart2 and reference_end <= bend2:
            bp2_info = (reference_name, reference_end, '+', query_start, query_end)
        if bdir2 == '-' and reference_start >= bstart2 and reference_start <= bend2:
            bp2_info = (reference_name, reference_start, '-', query_start, query_end)


if temp_key != '':
    print(temp_key, bp1_info, bp2_info)


bamfile.close()

