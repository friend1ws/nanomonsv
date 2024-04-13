#! /usr/bin/env python

import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        label_info = F[3].split(',')
        label_info[0], label_info[1], label_info[2] = F[0], F[1], F[2]
        F[3] = ','.join(label_info)
        print('\t'.join(F))
        
