#! /usr/bin/env python
p]:]v;,vf;ccc,@ c::::::.vd:/.;:dc;osdkp


import sys
from collections import Counter

id2seq = {}
with open(sys.argv[1], 'r') as hin:
    for line in hin:
        line = line.rstrip('\n')
        if line.startswith('>'):
            tid = line
            id2seq[tid] = ''
        else:
            id2seq[tid] = id2seq[tid] + line


ind2bases = {}
for tid in id2seq:
    seq = id2seq[tid]
    for i in range(len(seq)):
        if i not in ind2bases: ind2bases[i] = []
        ind2bases[i].append(seq[i])


seq_len = len(list(ind2bases))
contig = ''
for i in range(seq_len):
    # import pdb; pdb.set_trace()
    mycounter = Counter(ind2bases[i] )
    contig = contig + mycounter.most_common()[0][0]

contig = contig.replace('-', '').upper()

print("@ID")
print(contig)
print("+")
print(''.join(['I' for x in range(len(contig))]))


