#! /usr/bin/env python

import sys, gzip
import pysam

rmsk_file = sys.argv[1]

with gzip.open(rmsk_file, 'rt') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        # if len(F[5]) > 5: continue
        if int(F[2]) - int(F[1]) < 5800: continue
        if F[7] != "L1": continue
        if not F[3] in ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5"]: continue
        label = ','.join([F[0], F[1], F[2], F[5], F[3]])
        print('\t'.join([F[0], F[1], F[2], label, '0', F[5]]))
 
