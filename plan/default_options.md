# Make --use_racon and --single_bnd default on

## Status: TODO
- Created: 2026-03-09

## Overview
Change `--use_racon` and `--single_bnd` from opt-in to default-on in `nanomonsv get`.
Add `--no_racon` and `--no_single_bnd` flags for users who want the old behavior.

## Rationale
- `--use_racon` improves breakpoint accuracy and is recommended for all use cases.
- `--single_bnd` captures SVs at repeat boundaries that canonical SVs miss. Since v0.8.2, sbnd results have proper filtering (simple repeat, canonical SV overlap), making them reliable by default.
- Both options are already recommended in the tutorial. Making them default reduces user friction.

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
get.add_argument("--use_racon", default=True, action='store_true',
                 help="Use racon for error correction ... (default: True)")
get.add_argument("--no_racon", default=False, action='store_true',
                 help="Disable racon error correction")
get.add_argument("--single_bnd", default=True, action='store_true',
                 help="Generate single end breakpoints (default: True)")
get.add_argument("--no_single_bnd", default=False, action='store_true',
                 help="Disable single breakend SV detection")
```

### `nanomonsv/run.py`
Add logic to handle `--no_racon` / `--no_single_bnd`:
```python
if args.no_racon:
    args.use_racon = False
if args.no_single_bnd:
    args.single_bnd = False
```

### Documentation
- Update README usage examples (remove `--use_racon --single_bnd` from examples since they are now default)
- Update wiki Tutorial page
