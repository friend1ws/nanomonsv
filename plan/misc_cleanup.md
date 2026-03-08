# Clean up misc/ scripts

## Status: DONE
- Created: 2026-03-08

## Overview
Remove misc/ scripts whose functionality has been integrated into `nanomonsv get` as of v0.8.2.
Simplify remaining postprocess_sbnd pipeline.

## Removed files

### Simple repeat filtering (replaced by `nanomonsv get --simple_repeat_bed`)
- `misc/add_simple_repeat.py`
- `misc/subscript_postprocess_sbnd/add_simple_repeat.py`
- `misc/subscript_postprocess_sbnd/add_simple_repeat_sbnd.py`

### Post filtering (replaced by `Sv_filterer` in `post_proc.py`)
- `misc/post_filter_by_v0.3.py`

### Connect utilities (integrated into `nanomonsv/connect.py`)
- `misc/connect.py`
- `misc/get_edge.py`

## Simplified files

### postprocess_sbnd.sh
- Removed `add_simple_repeat.py` and `insert_classify` calls
- Now only runs: `annotate_contig.py` -> `plot_contig.R`

### integrate_sbnd.py -> annotate_contig.py
- Removed `Sv_filterer` / SV promotion logic (canonical SV overlap now handled by `nanomonsv get`)
- Removed `nanomonsv2info` canonical SV cross-reference
- Now generates: bwa.txt, rmsk.txt, class.txt

### Directory/file renames
- `misc/subscript_postprocess_sbnd/` -> `misc/subscript_sbnd/`
- `integrate_sbnd.py` -> `annotate_contig.py`
- `plot_sbnd_vis.R` -> `plot_contig.R`

## Kept files

### postprocess_sbnd pipeline
- `misc/postprocess_sbnd.sh` — Simplified pipeline script
- `misc/subscript_sbnd/annotate_contig.py` — BWA/RepeatMasker annotation + contig classification
- `misc/subscript_sbnd/plot_contig.R` — Visualization tool

### Utilities
- `misc/sv_type.py` — Adds SV_Type column (Del/Ins/Dup/Inv/Trans)
- `misc/make_test_bam.sh` — Test BAM generation (development use)

### Resources
- `misc/example/` — Sample result files (linked from README)
- `misc/img/` — Documentation images (used in wiki)
