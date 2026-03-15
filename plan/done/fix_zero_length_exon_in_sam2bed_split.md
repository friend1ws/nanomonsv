# Fix: ZeroDivisionError in insert_classify due to zero-length exon

## Status: DONE
- Completed: 2026-03-16
- Commit: 3693596

## Problem

`nanomonsv insert_classify` crashes with `ZeroDivisionError` at `insert_classify.py:247` (`pp_proc_filt_exon`).

```
ZeroDivisionError: float division by zero
  File "nanomonsv/insert_classify.py", line 247, in pp_proc_filt_exon
    match_ratio2 = float(F[12]) / (int(F[2]) - int(F[1]))
```

## Root Cause

`sam2bed_split()` generates zero-length exon records (start == end in BED) when minimap2 splice alignment produces a CIGAR with two consecutive N (intron) operators separated only by I (insertion), with no M (match) in between.

Example CIGAR from minimap2:
```
...58M1143N1628I25631N25M...
```

Processing flow in `sam2bed_split()`:
1. `58M` → `reference_pos_cur` advances by 58
2. `1143N` → exon is output, `reference_exon_cur` resets to `reference_pos_cur`
3. `1628I` → only `query_pos_cur` advances, `reference_pos_cur` does NOT change
4. `25631N` → exon is output, but `reference_pos_cur - reference_exon_cur = 0`

This produces a BED record like:
```
chr20	3855628	3855628	seq284,3693,551,2178	0	+
```

When this zero-length exon is later intersected with GENCODE exon annotations in `pp_proc_filt_exon()`, dividing by `(F[2] - F[1])` causes ZeroDivisionError.

## Actual Data

- Input: PAV insertion VCF (CHM13 reference), sample NA18939
- Problematic record: `chr1:54510521` INS 3693bp (SV_ID: `chr1-54510522-INS-3693`)
- minimap2 aligned this insertion sequence to chr20 with splice alignment mode (`-ax splice`)
- The CIGAR contained `...58M1143N1628I25631N...` pattern

minimap2 SAM (primary alignment):
```
seq284,3693  0  chr20  3854203  1  270S109M1D49M4D62M2I58M1143N1628I25631N25M1D58M4I2M1I2M2I138M8D8M3D2M905N24M2I32M4I3M8I3M1I2M1I4M1I2M4I5M3D2M3I2M4I4M4I11M9353N...
```

Resulting zero-length BED records:
```
chr20	3855628	3855628	seq284,3693,551,2178	0	+
```

## Fix

In `sam2bed_split()`, skip output when `exon_size == 0`:

```python
# Before
elif cigar[0] == 3:
    exon_size = reference_pos_cur - reference_exon_cur
    query_name2 = query_name + ',' + ...
    print('\t'.join([...]), file = hout)

# After
elif cigar[0] == 3:
    exon_size = reference_pos_cur - reference_exon_cur
    if exon_size > 0:
        query_name2 = query_name + ',' + ...
        print('\t'.join([...]), file = hout)
```

A zero-length exon carries no meaningful alignment information, so skipping it is the correct behavior.
