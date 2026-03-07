# VCF output: Fix issues in canonical SV VCF conversion

## Status: DONE
- Created: 2026-03-08
- Completed: 2026-03-08

## Overview
Review and fix issues in the current VCF output (`nanomonsv/vcf_convert.py`)
for compliance with VCF 4.3 specification.

## Issues found

### Issue 1: BND ALT field — REF base on wrong side when Dir_1 is `-`

**Code** (`vcf_convert.py` L124-128):
```python
tbracket = ']' if F["Dir_2"] == '+' else '['
if F["Dir_1"] == '+':
    talt1 = f'{tref1}{tsvinsseq}{tbracket}{tchrom2}:{tpos2}{tbracket}'
else:
    talt1 = f'{tbracket}{tchrom2}:{tpos2}{tbracket}{tsvinsseq}{tref2}'  # <-- uses tref2 instead of tref1
```

When Dir_1 is `-`, the ALT should use `tref1` (the REF base at this record's position),
not `tref2` (the REF base at the mate's position).

**Example from output** (r_0):
```
chr1  108952236  r_0_0  T  ]chr3:110694545]GCTTACCAGGT  ...
```
Here Dir_1=`-`, and the trailing base is `T`. This happens to match `tref1` (`T` at chr1:108952236),
but only by coincidence — the code actually uses `tref2`.

Per VCF spec: the reference base in ALT must be the base at POS of **this** record.

### Issue 2: BND mate record — same problem in reverse

**Code** (`vcf_convert.py` L139-144):
```python
tbracket = ']' if F["Dir_1"] == '+' else '['
if F["Dir_2"] == '+':
    talt2 = f'{tref2}{tsvinsseq}{tbracket}{tchrom1}:{tpos1}{tbracket}'
else:
    talt2 = f'{tbracket}{tchrom1}:{tpos1}{tbracket}{tsvinsseq}{tref2}'  # correct here (uses tref2)
```

The Dir_2=`-` case uses `tref2` which is correct for the mate record.
The Dir_2=`+` case uses `tref2` which is also correct.
So the mate record side appears correct. Only Issue 1 (the first record when Dir_1=`-`) is wrong.

### Issue 3: SVINSSEQ in BND includes REF base

**Example from output**:
```
chr1  108952236  r_0_0  REF=T  ALT=]chr3:110694545]GCTTACCAGGT  SVINSSEQ=GCTTACCAGG
```

The ALT `]chr3:110694545]GCTTACCAGGT` has 11 characters after the bracket,
but SVINSLEN=10 and SVINSSEQ has 10 characters. The extra character is the REF base `T`.
This is actually correct per spec — the REF base is part of ALT notation but not the insertion.
No issue here.

### Issue 4: INS records — END should equal POS

**Code** (`vcf_convert.py` L79,93):
```python
tend = int(F["Pos_2"]) - 1
tinfo = f"END={tend};SVTYPE=INS;SVINSLEN={tsvinslen};SVINSSEQ={tsvinsseq}"
```

Per VCF spec, for insertions END should equal POS (since the insertion happens at a single point).
Currently, `tend = Pos_2 - 1` which may differ from `Pos_1` if there's a small deletion component.

**Example from output**:
```
chr11  132054373  i_4  A  <INS>  ...  END=132054373;SVTYPE=INS;SVINSLEN=6048;...
```
In this case END=POS, which is correct. But this depends on Pos_2 = Pos_1 + 1 in the input.

### Issue 5: VCF records not sorted by position

VCF spec recommends (and many tools require) records sorted by CHROM+POS.
Currently, records are output in input file order, and BND pairs are adjacent
(record 1 at chr1, record 2 at chr3), which breaks sort order.

### Issue 6: Missing `Simple_repeat` FILTER definition

When `--simple_repeat_bed` is used, `Is_Filter` can contain `Simple_repeat`,
but there is no `##FILTER=<ID=Simple_repeat,...>` in the header.

**Example from output**:
```
chr1  117029127  d_0  A  <DEL>  .  Simple_repeat  ...
```
The `Simple_repeat` filter value appears in the data but has no header definition.

## Priority

| Issue | Severity | Notes |
|-------|----------|-------|
| 1 | Medium | Wrong REF base in BND ALT; may affect downstream tools |
| 5 | Medium | Unsorted VCF breaks many tools (bcftools, etc.) |
| 6 | Low | Missing FILTER header; easy fix |
| 4 | Low | Usually correct in practice |

## Related files
- `nanomonsv/vcf_convert.py` — VCF conversion code
- `workflow_test/output/guppy_6.1.2/COLO829/COLO829.nanomonsv.result.vcf` — Example output
