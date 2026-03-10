# Usage: python proc_hgsvc3.py <MEI_Callset.csv.gz>
# Extract LINE/L1 records from HGSVC3 MEI callset and output BED format.
# Strand is taken from PALMER_INFO if available, else from L1ME-AID_INFO.

import sys, csv, gzip, re

mei_file = sys.argv[1]

with gzip.open(mei_file, 'rt') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['TE_Designation'] != 'LINE/L1':
            continue

        chrom = row['CHROM']
        pos = int(row['POS'])
        start, end = str(pos), str(pos + 1)
        rec_id = row['ID']

        # Strand from PALMER_INFO if available, else from L1ME-AID_INFO
        palmer = row['PALMER_INFO']
        lme = row['L1ME-AID_INFO']
        if palmer != 'NONE':
            strand_m = re.search(r'STRAND=([^;]+)', palmer)
            strand = strand_m.group(1) if strand_m else '*'
        else:
            ori_m = re.search(r'Orientation:([^;]+)', lme)
            ori = ori_m.group(1) if ori_m else '*'
            strand = '+' if ori == '+' else '-' if ori == '-' else '*'

        name = 'HGSVC3_' + rec_id
        label = ','.join([chrom, start, end, strand, name])
        print('\t'.join([chrom, start, end, label, '0', strand]))
