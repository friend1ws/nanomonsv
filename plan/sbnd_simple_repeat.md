# sbnd: Integrate simple repeat filtering into `nanomonsv get`

## Status: TODO
- Created: 2026-03-08

## Overview
Move the sbnd simple repeat filtering (currently a post-processing step in `misc/postprocess_sbnd.sh`)
into the `nanomonsv get` command, so it runs automatically when `--simple_repeat_bed` is specified.

## Background

### Current sbnd simple repeat filtering (post-processing)
`misc/subscript_postprocess_sbnd/add_simple_repeat_sbnd.py` runs as a separate post-processing step:
- For each sbnd record, fetches simple repeats overlapping `Pos_1 ± 50bp` via `pysam.TabixFile`
- If `Pos_1` falls within `[repeat_start - 50, repeat_end + 50]`, adds the repeat coordinates as `Is_Simple_Repeat` column
- Otherwise, sets `Is_Simple_Repeat` to `"PASS"`

### Canonical SV simple repeat filtering (already integrated)
For canonical SVs (non-sbnd), `nanomonsv get` already performs simple repeat filtering via
`Sv_filterer.filter_simple_repeat()` in `nanomonsv/post_proc.py` (L177-195):
- Uses `--simple_repeat_bed` option
- Checks if **both** pos1 and pos2 fall within a simple repeat region (margin = 30bp)
- Appends `"Simple_repeat"` to the `Is_Filter` column (filter list approach)

### Key differences between canonical SV and sbnd
| | Canonical SV | sbnd (current script) |
|---|---|---|
| Breakpoints | Two (pos1, pos2) | One (pos1 only) |
| Logic | Both pos1 and pos2 must be within a repeat region | pos1 must be within a repeat region |
| Margin | 30bp | 50bp |
| Output | Appends to `Is_Filter` column | Adds `Is_Simple_Repeat` column with coordinates |

## Design

### Approach
Add simple repeat filtering inside `integrate_realignment_result_sbnd()` in `post_proc.py`.
When `--simple_repeat_bed` is provided, filter sbnd records using the same filter-list approach
as canonical SVs (appending `"Simple_repeat"` to `Is_Filter`).

### Changes required

#### 1. `nanomonsv/post_proc.py`
- Add a helper function `filter_sbnd_in_simple_repeat(tchr, tpos, simple_repeat_tb, margin=50)`:
  - Fetches simple repeats overlapping `pos ± margin` via `pysam.TabixFile.fetch()`
  - Returns `True` if `pos` falls within `[repeat_start - margin, repeat_end + margin]`
- Add `simple_repeat_bed` parameter to `integrate_realignment_result_sbnd()`
- Open TabixFile if `simple_repeat_bed` is provided
- Add `Is_Filter` column to the output header
- For each sbnd record, call the helper function and set `Is_Filter` to `"Simple_repeat"` or `"PASS"`

#### 2. `nanomonsv/run.py`
- Pass `args.simple_repeat_bed` to the `integrate_realignment_result_sbnd()` call (L423)

### Design decisions
- **Margin = 50bp**: Keep the existing sbnd margin (50bp) rather than canonical SV's 30bp,
  since sbnd has only one breakpoint and a wider margin is appropriate
- **Output format**: Use `Is_Filter` column with filter-list approach (consistent with canonical SV output),
  instead of the current `Is_Simple_Repeat` column with coordinate strings

## Related files
- `misc/postprocess_sbnd.sh` — Current post-processing pipeline (L14: calls add_simple_repeat_sbnd.py)
- `misc/subscript_postprocess_sbnd/add_simple_repeat_sbnd.py` — Current sbnd simple repeat script
- `nanomonsv/post_proc.py` L177-195 — Canonical SV simple repeat filter (`Sv_filterer.filter_simple_repeat()`)
- `nanomonsv/post_proc.py` L343-403 — sbnd integration function (`integrate_realignment_result_sbnd()`)
- `nanomonsv/run.py` L423-425 — Call site for `integrate_realignment_result_sbnd()`
- `nanomonsv/arg_parser.py` L79-80 — `--simple_repeat_bed` argument definition
