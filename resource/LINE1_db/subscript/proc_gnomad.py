#! /usr/bin/env python

import sys, gzip
import pysam

mei_file = sys.argv[1]

vcf_file = pysam.VariantFile(mei_file)
for record in vcf_file.fetch():

    # lchr = 'chr' + record.contig   
    lchr = record.contig
    lstart, lend = str(record.pos - 1), str(record.pos)
    lid = record.id
    # lstrand = record.info["MEINFO"][3]
    lstrand = '*'
    label = ','.join([lchr, lstart, lend, lstrand, lid])
    laf = 0
    # laf = float(record.info["AF"][0])
    # if laf < 0.01: continue

    print('\t'.join([lchr, lstart, lend, label, str(laf), lstrand]))
    
