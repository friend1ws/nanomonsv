# Prepare HGSVC3 mobile element insertion database for insert_classify

## Status: TODO
- Created: 2026-03-08

## Overview
Add HGSVC3 LINE1 data to the LINE1 database and remove duplicates between all sources.

Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0/

## Current LINE1 db sources and issues

### Sources in `LINE1.chm13v2.0.bed` (19,081 records, no dedup)
| Source | Records | Type | Coordinates |
|--------|---------|------|-------------|
| rmsk (RepeatMasker) | 5,139 | Reference full-length L1 (L1HS, L1PA2-5) | Wide intervals (start << end) |
| gnomAD-SV v3 | 13,290 | Non-reference insertion | Point (start, start+1) |
| 1000genomes (umary) | 652 | Non-reference insertion | Point (start, start+1) |

### Duplicate problem
- gnomAD vs 1000genomes: 424 of 652 overlap within 100bp (65%)
- gnomAD vs rmsk: 9 overlaps (non-ref insertion at edge of reference full-length L1)
- 1000genomes vs rmsk: 20 overlaps (same pattern)

### HGSVC3 MEI data to add
From `MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz`:
- LINE/L1 records: 1,649 total
- Full-length (SV_Length >= 5800): 565 (L1HS: 540, L1PA2: 25)

Subfamily detail (PALMER SUBFAM): L1Ta(568), L1Ta1d(386), L1Ambig(177), L1Ta1nd(154), L1Ta0(97), L1PreTa(28), L1Pa2(1)

## BED format
```
chr  start  end  info_string  score  strand
```
Where `info_string` = `chr,start,end,strand,name`

## Processing steps

### 1. Extract HGSVC3 LINE1 from CSV and convert to BED
```bash
zcat MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz \
  | awk -F',' 'NR>1 && $76 == "LINE/L1"' \
  | python3 hgsvc3_to_bed.py > hgsvc3.line1.chm13v2.0.bed
```
Extract CHROM, POS, and from PALMER_INFO: STRAND, SUBFAM.
Output: `chr  pos  pos+1  chr,pos,pos+1,strand,HGSVC3_SUBFAM_ID  0  strand`

### 2. Merge with dedup: `merge_bed_with_dedup()`
Priority order: rmsk > HGSVC3 > 1000genomes > gnomAD.
rmsk has the most reliable coordinates (reference-based), so it takes highest priority.
HGSVC3 is next because it has subfamily annotation (SUBFAM) and strand information.

Input: list of BED files in priority order (highest first).
Logic: for each file, add records only if no existing record is within 100bp on the same chromosome.

```python
def merge_bed_with_dedup(bed_files, window=100):
    """Merge multiple BED files, dropping duplicates by proximity.

    Args:
        bed_files: list of BED file paths in priority order (highest first)
        window: records within this distance on the same chr are considered duplicates

    Returns:
        list of BED lines (sorted by chr, start)
    """
    # collected records: {chr: sorted list of (start, end, line)}
    records = {}

    for bed_file in bed_files:
        for line in open(bed_file):
            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            chr_records = records.setdefault(chrom, [])
            # binary search for nearby records
            if not _has_nearby(chr_records, start, end, window):
                bisect.insort(chr_records, (start, end, line.strip()))

    # flatten and sort
    result = []
    for chrom in sorted(records.keys()):
        for start, end, line in records[chrom]:
            result.append(line)
    return result


def _has_nearby(chr_records, start, end, window):
    """Check if any existing record overlaps [start-window, end+window]."""
    # binary search to find candidates
    idx = bisect.bisect_left(chr_records, (start - window,))
    for i in range(max(0, idx - 1), len(chr_records)):
        s, e, _ = chr_records[i]
        if s > end + window:
            break
        # two intervals overlap if NOT (e2+w < s1 or s2-w > e1)
        if not (e + window < start or s - window > end):
            return True
    return False
```

Usage:
```python
lines = merge_bed_with_dedup([
    "rmsk.line1.chm13v2.0.bed",
    "hgsvc3.line1.chm13v2.0.bed",
    "1000genomes.line1.chm13v2.0.bed",
    "gnomad.line1.chm13v2.0.bed",
])
```

### 3. Write, compress, and index
```bash
bgzip LINE1.chm13v2.0.bed
tabix -p bed LINE1.chm13v2.0.bed.gz
```

### Record counts (chm13v2.0)
| Source | Before dedup | After dedup |
|--------|-------------|-------------|
| rmsk | 5,139 | 5,134 |
| HGSVC3 | 1,664 | 1,646 |
| 1000genomes | 652 | 394 |
| gnomAD-SV v3 | 13,290 | 12,712 |
| **Total** | **20,745** | **19,886** |

### Final results (all references)
| Reference | Records |
|-----------|---------|
| hg38 | 20,006 |
| hg19 | 19,674 |
| chm13v2.0 | 19,886 |

## Same process for hg38 and hg19
- hg38: Apply the same steps using the GRCh38 callset (`MEI_Callset_GRCh38.ALL.20241211.csv.gz`).
- hg19: HGSVC3にはhg19データがないので、hg38のBEDからliftOverで変換する。

## Related files
- `nanomonsv/insert_classify.py` L529-658 — `organize_info()` uses LINE1_db
- `nanomonsv/data/LINE1.hg38.bed.gz` — Current bundled LINE1 db (8,490 records)
- `resource/LINE1_db/LINE1.hg38.bed` — Newer LINE1 db with gnomAD v3 (19,254 records)
- `resource/LINE1_db/MEI_Callset_GRCh38.ALL.20241211.csv.gz` — Downloaded HGSVC3 data
- `resource/LINE1_db/MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz` — Downloaded HGSVC3 data
