#! /usr/bin/env bash

import sys, csv
import pysam

input_file = sys.argv[1]
repeat_file = sys.argv[2]


repeat_tb = pysam.TabixFile(repeat_file)

with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n')
    print(f'{header}\tIs_Simple_Repeat')
    for line in hin:
        F = line.rstrip('\n').split('\t')
        
        # check simplerepeat annotation
        repeat = "PASS"
        tabixErrorFlag = 0
        try:
            records = repeat_tb.fetch(F[0], int(F[1]) - 50, int(F[1]) + 50)
        except Exception as inst:
            print("%s: %s" % (type(inst), inst.args), file = sys.stderr)
            tabixErrorFlag = 1

        if tabixErrorFlag == 0:
            for record_line in records:
                record = record_line.split('\t')
                if int(F[1]) >= int(record[1]) - 50 and int(F[1]) <= int(record[2]) + 50:
                    repeat = record[0] + ':' + record[1] + '-' + record[2]

        print('\t'.join(F) + '\t' + repeat)


