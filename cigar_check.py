#! /usr/bin/env python

import pysam

def parse_alignment_info(input_bam, output_file):

    hout = open(output_file, 'w') 
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch():

        # print(read)
        query_name = read.query_name
        query_start = str(read.query_alignment_start + 1)
        query_end = str(read.query_alignment_end)
        reference_name = read.reference_name
        reference_start = str(read.reference_start + 1)
        reference_end = str(read.reference_end)
        mapping_quality = str(read.mapping_quality)
        is_secondary = str(read.is_secondary)
        is_supplementary = str(read.is_supplementary)

        cigar_stats = read.get_cigar_stats()
        num_M = cigar_stats[0][0]
        num_I = cigar_stats[0][1]
        num_D = cigar_stats[0][2]

        num_M2, num_I2, num_D2 = 0, 0, 0
        cigartuples = read.cigartuples
        for elm in cigartuples:
            if elm[0] == 0: num_M2 = num_M2 + elm[1]
            if elm[0] == 1: num_I2 = num_I2 + elm[1]
            if elm[0] == 2: num_D2 = num_D2 + elm[1]
            
            if elm[0] == 1 and elm[1] >= 5000:
                print(reference_name + " Ins: " + str(elm[1]))
            if elm[0] == 2 and elm[1] >= 5000:
                print(reference_name + " Del: " + str(elm[1]))

        if num_M != num_M2 or num_I != num_I2 or num_D != num_D2:
            print("error")
    
        
    hout.close()
    bamfile.close()


if __name__ == "__main__":

    import sys

    input_bam = sys.argv[1]
    output_file = sys.argv[2]
    parse_alignment_info(input_bam, output_file)

