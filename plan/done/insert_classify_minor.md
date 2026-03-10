# insert_classify: minor fixes

## Status: DONE
- Created: 2026-03-09
- Completed: 2026-03-10 (commit 2e9bff0)

## Overview
Minor code quality fixes in `nanomonsv/insert_classify.py`.

## Issues

### 1. Duplicate file opens (resource leak)

#### `summarize_bwa_alignment` (L438, L448)
`samfile = pysam.AlignmentFile(input_sam, 'r')` is opened twice. L438 is unused.
- Fix: Remove L438.

#### `summarize_bwa_alignment2` (L488, L499)
Same issue. `samfile` opened twice, L488 is unused.
- Fix: Remove L488.

#### `summarize_bwa_alignment2` (L505, L571)
`hout = open(output_file, 'w')` opened at L505, never used or closed. Reopened at L571.
- Fix: Remove L505.

### 2. Unused variables in `sam2bed_split` (L140-141)
`is_secondary` and `is_supplementary` are assigned but never used. The code uses `read.is_secondary` and `read.is_supplementary` directly.
- Fix: Remove L140-141.

### 3. Typos in `annotate_sv_file` header (L753)
- `"L1_Raito"` → `"L1_Ratio"`
- `"SV_Ratio"` → `"SVA_Ratio"`
- `"PDS_Exon_Num"` → `"PSD_Exon_Num"`
- Fix: Correct the header strings.

## Related files
- `nanomonsv/insert_classify.py`
