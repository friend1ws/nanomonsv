#! /usr/bin/env python

import sys, gzip
import pysam

rmsk_file = sys.argv[1]

with gzip.open(rmsk_file, 'rt') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        # if len(F[5]) > 5: continue
        if int(F[7]) - int(F[6]) < 5800: continue
        if F[12] != "L1": continue
        label = ','.join([F[5], F[6], F[7], F[9], F[10]])
        print('\t'.join([F[5], F[6], F[7], label, '0', F[9]]))
 
