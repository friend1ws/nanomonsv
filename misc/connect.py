#! /usr/bin/env python

import sys

input_file = sys.argv[1]

segments = []
with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        segments.append([F[:16] + ['+']])



update = True
while update:
    update = False
    for seg1 in range(len(segments) - 1):
        for seg2 in range(seg1 + 1, len(segments)):

            break1_1 = '\t'.join(segments[seg1][0][:3])
            break1_2 = '\t'.join(segments[seg1][0][4:7])
            break1_3 = '\t'.join(segments[seg1][-1][7:10])
            break1_4 = '\t'.join(segments[seg1][-1][11:14])

            break2_1 = '\t'.join(segments[seg2][0][:3])
            break2_2 = '\t'.join(segments[seg2][0][4:7])
            break2_3 = '\t'.join(segments[seg2][-1][7:10])
            break2_4 = '\t'.join(segments[seg2][-1][11:14])

            if break1_3 == break2_1 and break1_4 == break2_2:
                sdir1, sdir2 = '+', '-'
                update = True
            elif break1_3 == break2_4 and break1_4 == break2_3:
                sdir1, sdir2 = '+', '+'
                update = True
            elif break1_1 == break2_3 and break1_2 == break2_4:
                sdir1, sdir2 = '-', '+'
                update = True
            elif break1_1 == break2_2 and break1_2 == break2_1:
                sdir1, sdir2 = '-', '-'
                update = True

            if update: break

        if update: break


    if update:
        if sdir1 == sdir2:
            for j in range(len(segments[seg2])):
                segments[seg2][j] = segments[seg2][j][11:14] + [segments[seg2][j][10]] + segments[seg2][j][7:10] + \
                                    segments[seg2][j][4:7] + [segments[seg2][j][3]] + segments[seg2][j][:3] + \
                                    segments[seg2][j][14:16] + ['+' if segments[seg2][j][16] == '-' else '-']

                segments[seg2].reverse()

        if sdir1 == '+':
            segments[seg1] = segments[seg1] + segments[seg2]
        else:
            segments[seg1] = segments[seg2] + segments[seg1]

        del(segments[seg2])
        

for seg in segments:

    for s in seg:

        print(s)


    print()
    


"""
8       130933956       +       ---     8       130871365       -       8       130872760       +       ---     8       131401249       +       True    1395 
"""



