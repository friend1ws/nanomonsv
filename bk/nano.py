#! /usr/bin/env python

import pysam

def parse_alignment_info(input_bam, output_file):

    hout = open(output_file, 'w') 
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch():

        # print(read)
        if read.is_secondary: continue

        query_name = read.query_name
        query_strand = '-' if read.is_reverse else '+'
        query_length = read.infer_read_length()

        reference_name = read.reference_name
        reference_start = str(read.reference_start + 1)
        reference_end = str(read.reference_end)
        mapping_quality = str(read.mapping_quality)
        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary
        
        cigar_stats = read.get_cigar_stats()
        num_M = cigar_stats[0][0]
        num_I = cigar_stats[0][1]
        num_D = cigar_stats[0][2]

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
            
        print('\t'.join([query_name, str(query_start), str(query_end), str(query_length), query_strand, \
                         reference_name, reference_start, reference_end, mapping_quality, \
                         str(num_M), str(num_I), str(num_D), str(is_supplementary), str(is_secondary)]), file = hout)

    hout.close()
    bamfile.close()


if __name__ == "__main__":

    import sys

    input_bam = sys.argv[1]
    output_file = sys.argv[2]
    parse_alignment_info(input_bam, output_file)

