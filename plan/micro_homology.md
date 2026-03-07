# Micro homology detection and reporting

## Status: TODO
- Created: 2026-03-08

## Overview
Detect micro homology at SV breakpoints and report it in the output,
particularly in VCF INFO fields (HOMLEN, HOMSEQ).

## Background

### What is micro homology at breakpoints?
When a structural variant has identical sequences at both sides of the breakpoint,
the exact breakpoint position is ambiguous. For example:

```
Reference region1: ...ATCGATCG|NNNNN...    (breakpoint after G)
Reference region2: ...NNNNN|ATCGATCG...    (breakpoint before A)
Contig:            ...ATCGATCGATCG...
```

If the contig matches `ATCG` at both sides, the breakpoint could be placed
at any position within those 4 bases. This shared sequence is "micro homology."

### VCF INFO fields (VCF 4.3 spec)
- `HOMLEN` — Length of micro homology at breakpoint (Integer)
- `HOMSEQ` — Micro homology sequence (String)
- `CIPOS` — Confidence interval around POS (Integer pair, e.g., `CIPOS=-4,0` for left-aligned)

Example VCF record with micro homology:
```
chr1  12345  sv_0  A  <DEL>  .  PASS  SVTYPE=DEL;END=56789;SVLEN=-44444;HOMLEN=4;HOMSEQ=ATCG;CIPOS=0,4
```

### When breakpoint is left-aligned:
- POS is at the leftmost valid position
- CIPOS=0,HOMLEN (uncertainty extends rightward)
- This is the conventional choice for VCF (consistent with variant left-alignment convention)

## Current implementation

### smith_waterman.py — `sw_jump()`
A custom two-region Smith-Waterman with a jump from H1 (region1) to H2 (region2).

Key return values: `(score, contig_align, region1_align, region2_align, contig_seq, region_seq)`
- `contig_align = (i1_start, i1_end, i2_start, i2_end)` — contig coordinates for each region
- The relationship between `i1_end` and `i2_start` determines micro homology:
  - `i1_end < i2_start`: inserted sequence exists between the two alignments
  - `i1_end + 1 == i2_start`: clean breakpoint (no micro homology, no insertion)
  - `i1_end >= i2_start`: **micro homology** of length `i1_end - i2_start + 1`

### locate_bp.py — `get_refined_bp()`
Calls `sw_jump()` and converts alignment coordinates to genomic breakpoint positions.

**Current behavior** (L44-53):
```python
bp_pos1 = bstart1 + region1_align[1] - 1 if dir1 == '+' else bend1 - region1_align[1] + 1
bp_pos2 = bstart2 + region2_align[0] - 1 if dir2 == '-' else bend2 - region2_align[0] + 1

if contig_align[2] - contig_align[1] == 1:
    inseq = '---'
elif contig_align[2] - contig_align[1] > 1:
    inseq = contig[(contig_align[1]):(contig_align[2] - 1)]
```

- Uses `region1_align[1]` (end of region1 alignment) and `region2_align[0]` (start of region2 alignment)
  to determine one fixed breakpoint position
- When `contig_align[1] >= contig_align[2]` (i.e., overlap/micro homology), the code falls through
  to the `else` branch which only logs a warning ("Alignment inconsistent!!")
- **Micro homology information is lost** — not stored or reported

### Problem
1. The traceback in `sw_jump()` finds a single optimal path, choosing one jump point
2. When micro homology exists, there are multiple equally valid jump points
   (each corresponding to a different left/right-aligned breakpoint)
3. Currently, only one breakpoint is reported and micro homology is not detected

## Design

### Approach

#### Step 1: Detect micro homology from sw_jump result
After `sw_jump()` returns, compare the contig sequence around the jump point
with both reference regions to determine the extent of micro homology.

Specifically, given:
- `i1_end` — last contig position aligned to region1
- `i2_start` — first contig position aligned to region2

If `i1_end >= i2_start`, the overlapping contig bases are the micro homology sequence.

Even when `i1_end < i2_start` (insertion case), there may still be micro homology
if the bases adjacent to the breakpoint match on both sides. Need to extend the comparison
to check for this.

#### Step 2: Report micro homology length and sequence
Add to the refined breakpoint output:
- Micro homology length
- Micro homology sequence

#### Step 3: Left-align breakpoints
Choose the leftmost valid breakpoint position (consistent with VCF convention).
Report `CIPOS=0,HOMLEN` to indicate the ambiguity range.

#### Step 4: Add VCF INFO fields
In `vcf_convert.py`, output:
- `HOMLEN=<length>` — micro homology length (0 if none)
- `HOMSEQ=<sequence>` — micro homology sequence (omitted if none)
- `CIPOS=0,<HOMLEN>` — confidence interval for POS

### Open questions
- Should the existing `sw_jump()` be modified to enumerate all optimal paths,
  or is post-hoc micro homology detection (comparing sequences around the jump point) sufficient?
- How should micro homology interact with inserted sequence?
  (e.g., complex breakpoints with both micro homology and insertion)
- Performance impact of enumerating paths vs post-hoc detection?

## Related files
- `nanomonsv/smith_waterman.py` — SW alignment with jump (`sw_jump()`)
- `nanomonsv/locate_bp.py` L13-59 — `get_refined_bp()` converts alignment to breakpoint positions
- `nanomonsv/locate_bp.py` L89-109 — `locate_bp()` main loop
- `nanomonsv/vcf_convert.py` — VCF output (needs HOMLEN/HOMSEQ/CIPOS)
