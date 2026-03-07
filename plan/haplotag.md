# Haplotype-aware read counting and filtering

## Status: TODO
- Created: 2026-03-08

## Overview
When BAM reads have haplotype tags (HP:i:1, HP:i:2), count supporting reads
per haplotype and filter SVs where supporting reads are dispersed across
both haplotypes (likely artifacts).

## Background

### HP tag in BAM
Long-read phasing tools (e.g., WhatsHap, LongPhot, margin) add the `HP` tag:
- `HP:i:1` — assigned to haplotype 1
- `HP:i:2` — assigned to haplotype 2
- No HP tag — unphased

A true somatic SV should have supporting reads predominantly from one haplotype
(or from unphased reads). If supporting reads are split roughly evenly between
H1 and H2, the call is likely an artifact (e.g., alignment error at a
heterozygous germline variant).

## Current implementation

### `gather_local_read_for_realignment()` (L27-134)
Fetches reads from BAM around breakpoint positions. For each read, stores:
- `read.qname` — read name
- `read.mapping_quality` — mapping quality
- `read.query_sequence` — read sequence

**HP tag is not read.** The read object has `read.get_tag("HP")` available
but it is never called.

### Intermediate file format
```
key\treadid\tmapq1\tmapq2\tseq
```
No haplotype information is included.

### `Alignment_counter.count_alignment_and_print()` (L292-358)
Counts total `all_rnames` and `supporting_reads`. Outputs:
```
chr1\tpos1\tdir1\tchr2\tpos2\tdir2\tinseq\tid\tall_count\tsupporting_count
```
No per-haplotype breakdown.

### `Alignment_counter.count_alignment_and_print_sbnd()` (L361-397)
Same pattern for single breakends. No haplotype awareness.

## Design

### Approach
Always collect HP tag from BAM reads. Report per-haplotype supporting read
counts. When HP tags are present, add a filter for haplotype dispersion.

### Changes required

#### 1. `gather_local_read_for_realignment()` — Extract HP tag
When iterating reads from BAM, extract HP tag:
```python
try:
    hp = read.get_tag("HP")
except KeyError:
    hp = 0  # unphased
```

Store HP alongside mapq in `key2rname2mapq`:
- Current: `key2rname2mapq[key][read.qname] = [mapq1, mapq2]`
- New: `key2rname2mapq[key][read.qname] = [mapq1, mapq2, hp]`

For sbnd: `key2rname2mapq_sbnd[key][read.qname] = [mapq, hp]`

Note: A read is fetched from two regions (around pos1 and pos2). The HP tag
should be the same regardless of which fetch returns it, so storing it from
either fetch is fine.

#### 2. Intermediate file format
Add HP as a new column:
```
key\treadid\tmapq1\tmapq2\thp\tseq
```

#### 3. `Alignment_counter` — Store per-read haplotype
Add `self.readid2hp = {}` alongside existing `self.readid2mapq`.
Modify `add_long_read_seq()` to accept and store HP:
```python
def add_long_read_seq(self, treadid, tseq, tmapq1, tmapq2, thp):
    self.readid2mapq[treadid] = [...]
    self.readid2hp[treadid] = int(thp)
```

#### 4. `count_alignment_and_print()` — Per-haplotype counts
After determining `supporting_reads`, split by haplotype:
```python
hp1_reads = [r for r in supporting_reads if self.readid2hp[r] == 1]
hp2_reads = [r for r in supporting_reads if self.readid2hp[r] == 2]
hp0_reads = [r for r in supporting_reads if self.readid2hp[r] == 0]
```

Output additional columns:
```
...\tall_count\tsupporting_count\thp1_count\thp2_count\thp0_count
```

Same for `count_alignment_and_print_sbnd()`.

#### 5. Haplotype dispersion filter
In `post_proc.py`, after reading supporting read counts:
- If HP tags are present (hp1_count + hp2_count > 0):
  - Calculate dominant haplotype ratio: `max(hp1, hp2) / (hp1 + hp2)`
  - If ratio < threshold (e.g., 0.9), flag as `Haplotype_Dispersion`
- If no HP tags (hp1_count + hp2_count == 0), skip the filter

This is a filter flag (like Simple_repeat), not an outright rejection.

#### 6. Output format
Add columns to result.txt:
- `Hp1_Support_Num` — supporting reads from haplotype 1
- `Hp2_Support_Num` — supporting reads from haplotype 2
- `Hp0_Support_Num` — supporting reads unphased

Add `Haplotype_Dispersion` to Is_Filter when applicable.

### Open questions
- Should the threshold be configurable via command-line argument?
  (e.g., `--hp_ratio_thres 0.9`)
- How to handle cases where most supporting reads are unphased (hp0)?
  Should unphased reads be excluded from the ratio calculation?
- Should control BAM also be checked for haplotype consistency?
- Should per-haplotype counts be reported in VCF INFO fields?
  (e.g., `HP1_SR=3;HP2_SR=1;HP0_SR=2`)

## Related files
- `nanomonsv/count_sread_by_alignment.py` L27-134 — `gather_local_read_for_realignment()`
- `nanomonsv/count_sread_by_alignment.py` L167-398 — `Alignment_counter` class
- `nanomonsv/count_sread_by_alignment.py` L292-358 — `count_alignment_and_print()`
- `nanomonsv/count_sread_by_alignment.py` L361-397 — `count_alignment_and_print_sbnd()`
- `nanomonsv/post_proc.py` — Where filtering decisions are made
