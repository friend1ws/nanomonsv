# Make racon and --single_bnd default on

## Status: DONE
- Created: 2026-03-09
- Completed: 2026-03-16
- Note: internal variable name `use_racon` remains in generate_consensus.py (cosmetic only; logic is correct via `args.use_racon = not args.use_mafft`)

## Overview
Change `--use_racon` and `--single_bnd` from opt-in to default-on in `nanomonsv get`.
- Racon becomes the default consensus method. Add `--use_mafft` for users who want the old mafft-based behavior.
- `--single_bnd` becomes default-on. Add `--no_single_bnd` to disable.
- Remove `--use_racon` flag (no longer needed since racon is default).

## Rationale
- Racon improves breakpoint accuracy and is recommended for all use cases.
- `--single_bnd` captures SVs at repeat boundaries that canonical SVs miss. Since v0.8.2, sbnd results have proper filtering (simple repeat, canonical SV overlap), making them reliable by default.
- Both options are already recommended in the tutorial. Making them default reduces user friction.
- `--use_mafft` is clearer than `--no_racon` because it tells the user what will actually be used.

## Changes required

### `nanomonsv/arg_parser.py` (L139-143)

Current:
```python
get.add_argument("--use_racon", default=False, action='store_true',
                 help="Use racon for error correction ... (default: False)")
get.add_argument("--single_bnd", default=False, action='store_true',
                 help="Generate single end breakpoints (default: False)")
```

New:
```python
get.add_argument("--use_mafft", default=False, action='store_true',
                 help="Use mafft instead of racon for consensus generation")
get.add_argument("--single_bnd", default=True, action='store_true',
                 help="Generate single end breakpoints (default: True)")
get.add_argument("--no_single_bnd", default=False, action='store_true',
                 help="Disable single breakend SV detection")
```

### `nanomonsv/run.py`
- Replace `args.use_racon` with `not args.use_mafft` (racon is now default, mafft is opt-in)
- Handle `--no_single_bnd`:
```python
if args.no_single_bnd:
    args.single_bnd = False
```
- Update tool check logic:
```python
if args.use_mafft:
    is_tool("mafft")
else:
    is_tool("racon")
```
- Remove the `single_bnd` requires `use_racon` check (racon is now always used unless `--use_mafft`)

### `nanomonsv/generate_consensus.py`
- Rename `use_racon` parameter to `use_mafft` (invert logic) throughout

### Documentation
- Update README usage examples (remove `--use_racon --single_bnd` from examples since they are now default)
- Update wiki Tutorial page
