#! /usr/bin/env python

import pysam

def parse_alignment_info(input_bam, output_file):

    hout = open(output_file, 'w') 
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    for read in bamfile.fetch():

        if read.is_secondary: continue
        if read.is_supplementary: continue

        cigartuples = read.cigartuples
        Is = [x[1] for x in cigartuples if x[0] == 1]
        Ds = [x[1] for x in cigartuples if x[0] == 2]
        if (len(Is) == 0 or max(Is) <= 500) and (len(Ds) == 0 or max(Ds) <= 500): continue

        print(read.query_name, read.reference_name + '\t' + str(read.reference_start) + '\t' + str(max(Is)) + '\t' + str(max(Ds)))

        """"
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

        print('\t'.join([query_name, query_start, query_end, reference_name, reference_start, reference_end, mapping_quality, \
                        str(num_M), str(num_I), str(num_D), is_supplementary, is_secondary]), file = hout)
        """
    
    hout.close()
    bamfile.close()


if __name__ == "__main__":

    import sys

    input_bam = sys.argv[1]
    output_file = sys.argv[2]
    parse_alignment_info(input_bam, output_file)

#! /usr/bin/env python


