#! /usr/bin/env python

import pysam

bamfile = pysam.AlignmentFile("../../workspace/seq_data/nanopore/HCC1954.bam", 'rb')

# for read in bamfile.fetch("12", 41451383, 41451403):
#     if read.query_name == "00001773-fff3-4751-b8db-5339de63362c":
for read in bamfile.fetch("22", 36372090, 36372110):
    if read.query_name == "0000232e-568e-4fc8-8c22-3a4aeb13ff7b":
        import pdb; pdb.set_trace()
        print(read)

