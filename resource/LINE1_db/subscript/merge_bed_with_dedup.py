# Usage: python merge_bed_with_dedup.py <bed1> <bed2> [bed3] [bed4] ...
# Merge multiple BED files with priority-based dedup.
# Files are given in priority order (highest first).
# A record is dropped if any higher-priority record is within 100bp on the same chr.

import sys, bisect

def merge_bed_with_dedup(bed_files, window=100):
    # {chr: sorted list of (start, end, line)}
    records = {}

    for bed_file in bed_files:
        for line in open(bed_file):
            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            chr_records = records.setdefault(chrom, [])
            if not _has_nearby(chr_records, start, end, window):
                bisect.insort(chr_records, (start, end, line.strip()))

    for chrom in sorted(records.keys()):
        for start, end, line in records[chrom]:
            print(line)


def _has_nearby(chr_records, start, end, window):
    idx = bisect.bisect_left(chr_records, (start - window,))
    for i in range(max(0, idx - 1), len(chr_records)):
        s, e, _ = chr_records[i]
        if s > end + window:
            break
        if not (e + window < start or s - window > end):
            return True
    return False


merge_bed_with_dedup(sys.argv[1:])
