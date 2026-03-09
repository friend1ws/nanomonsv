# sbnd: Remove single breakends overlapping with canonical SVs

## Status: DONE
- Created: 2026-03-08
- Completed: 2026-03-08

## Overview
Remove single breakend (sbnd) calls whose breakpoint overlaps with a canonical SV breakpoint.
This is currently done in the post-processing script `misc/subscript_postprocess_sbnd/integrate_sbnd.py`,
and should be integrated into `nanomonsv get`.

## Analysis of integrate_sbnd.py

`integrate_sbnd.py` is a large post-processing script that does multiple things.
Below is a breakdown of each step.

### Step 1: Build canonical SV breakpoint list (L258-265)
Reads `{prefix}.nanomonsv.result.txt` (canonical SV results) and builds a dict `nanomonsv2info`
keyed by `(Chr, Pos, Dir)` for both breakpoints of each SV.

### Step 2: Filter sbnd by canonical SV overlap + generate contig FASTA (L267-271)
Reads `{prefix}.nanomonsv.sbnd.result.txt` and for each sbnd record:
- Calls `key_check()` to check if `(Chr_1, Pos_1, Dir_1)` matches any canonical SV breakpoint
  within margin=50bp
- If **no overlap**: writes the contig to a FASTA file for downstream analysis
- If **overlap**: skips it (effectively removing sbnd entries that duplicate canonical SVs)

### Step 3: RepeatMasker + BWA alignment of contigs (L273-286)
- Runs RepeatMasker on the filtered contigs
- Runs `bwa mem` to align contigs to the reference genome
- Parses RepeatMasker output (`proc_rmsk`) and SAM output (`proc_sam`)

### Step 4: Classify contigs (L288-334)
For each contig, creates a `Contig_info` object and classifies it:
- `Simple_repeat` — >80% of contig covered by simple repeats
- `Satellite` / `Satellite/centr` — >80% covered by satellite sequences
- `Plain_SV` — contig has a BWA alignment >=2000bp with MQ>=40 starting within first 100bp
- `L1_Mediated_Del` — L1HS/L1P1/L1PA2 segments followed by a long BWA alignment
- `Complex` — none of the above

### Step 5: Re-process canonical SVs + add sbnd-derived SVs (L337-360)
- Creates a new `Sv_filterer` and adds all canonical SVs
- For sbnd contigs classified as `Plain_SV` or `L1_Mediated_Del`, extracts SV coordinates
  from the contig alignment and adds them to the filterer as additional SVs
- Runs `apply_filters()` and writes `{prefix}.nanomonsv.proc.result.txt`
  (a merged, re-filtered canonical SV result that includes sbnd-derived SVs)

### Step 6: Write filtered sbnd output (L363-374)
- Reads sbnd result again, excludes contigs classified as `Plain_SV` or `L1_Mediated_Del`
  (since those were promoted to canonical SVs in Step 5)
- Adds `SBND_Class` column and writes `{prefix}.nanomonsv.sbnd.proc.result.txt`

## What integrate_realignment_result_sbnd() currently does (post_proc.py L343-403)

The function in `nanomonsv get` already:
1. **Builds `nanomonsv_bp_list`** (L347-354): reads canonical SV result file and collects
   breakpoint regions `(chr, pos-30, pos+30)` for both breakpoints of each SV.
2. **But the overlap check is commented out** (L379-384): the code that skips sbnd records
   overlapping with canonical SV breakpoints is wrapped in a triple-quoted string,
   so it is currently **inactive**.
3. The rest of the function filters sbnd by tumor/control read counts and VAF, then writes output.

### Differences between integrate_sbnd.py and current get
| Aspect | integrate_sbnd.py | current get (post_proc.py) |
|--------|-------------------|---------------------------|
| Overlap check | Active, margin=50bp, checks (chr, pos, dir) | Commented out, margin=30bp, checks only (chr, pos) without dir |
| What happens to overlapping sbnd | Excluded from output | Nothing (check disabled) |
| Contig classification | Yes (RepeatMasker + BWA) | No |
| Plain_SV/L1 promotion to canonical | Yes | No |
| SBND_Class column | Yes | No |

## Design

### Scope of this task
Only the **overlap removal** part (Steps 1-2 of integrate_sbnd.py).
Contig classification, RepeatMasker/BWA analysis, and SV promotion (Steps 3-6) are out of scope.

The overlap removal logic is purely coordinate-based (`key_check()` compares chr/pos/dir within a margin).
It has **no dependency on RepeatMasker or BWA** — those are only used in Steps 3-6 for contig classification.

### Approach
Uncomment and fix the existing overlap check code in `integrate_realignment_result_sbnd()`:

1. The `nanomonsv_bp_list` construction (L347-354) already exists and works.
   - Consider whether to also store `Dir` for the check (integrate_sbnd.py does, current code does not).
   - Consider margin: integrate_sbnd.py uses 50bp, current code uses 30bp.

2. Uncomment L379-384 and activate the overlap check.

### Open questions
- **Should Dir be checked?** integrate_sbnd.py checks (chr, pos, dir) match,
  current commented-out code only checks (chr, pos). Checking dir is more precise.
- **Margin value**: 50bp (integrate_sbnd.py) vs 30bp (current code). Should align with one.
- **Output format**: Should filtered sbnd records be excluded entirely (like integrate_sbnd.py)
  or marked in an `Is_Filter` column (like canonical SV filtering)?

## Related files
- `misc/subscript_postprocess_sbnd/integrate_sbnd.py` — Current post-processing script (full pipeline)
- `misc/postprocess_sbnd.sh` L8 — Calls integrate_sbnd.py
- `nanomonsv/post_proc.py` L343-403 — `integrate_realignment_result_sbnd()` with commented-out overlap check
- `nanomonsv/run.py` L423-425 — Call site for `integrate_realignment_result_sbnd()`
