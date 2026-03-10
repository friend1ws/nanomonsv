# insert_classify: annotation logic audit

## Status: TODO
- Created: 2026-03-10

## Overview
Audit of the annotation logic in `nanomonsv/insert_classify.py`, covering data flow, classification decisions, and identified issues.

## Data flow

```
inserted sequence
    ├── RepeatMasker  → summarize_rmsk()           → repeat_type, L1/Alu/SVA_ratio, rmsk_info
    ├── BWA alignment → summarize_bwa_alignment2()  → alignment_infos, inserted_pos
    └── check_tsd_polyAT()                          → is_polyAT, tsd
            │
            ▼
    organize_info()  ← LINE1_db
        → source_info[0..10]:
            [0] repeat_type        "Simple_LINE1", "Inverted_LINE1", "Other_LINE1", "Alu", "SVA", "None"
            [1] l1_ratio           float as str
            [2] alu_ratio          float as str
            [3] sva_ratio          float as str
            [4] rmsk_info          semicolon-separated repeat entries (from summarize_rmsk)
            [5] alignment_infos    semicolon-separated alignment entries
            [6] inserted_pos       "chr,start,end" or "---"
            [7] is_polyAT          "polyA", "polyT", "None", or "---"
            [8] tsd                sequence or "None" or "---"
            [9] line1_info         LINE1 source ID or "None" (Python str(None))
           [10] transduction_class "Solo", "Partnered", "Orphan", or "None" (Python str(None))
            │
            ▼
    annotate_sv_file() / annotate_vcf_file()
        → insert_type:   Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD, "---"
        → is_inversion:  Simple, Inverted, Other, NA
        + ppseudo_info:  PSD_Gene, PSD_Overlap_Ratio, PSD_Exon_Num
```

## Classification logic

### repeat_type (proc_rmsk_info, L317-337)
- L1_ratio >= 0.8: 1 L1 segment → "Simple_LINE1", 2 segments opposite strand → "Inverted_LINE1", else → "Other_LINE1"
- Alu_ratio >= 0.8 → "Alu"
- SVA_ratio >= 0.8 → "SVA"
- Otherwise → "None"
- Note: Alu/SVA checks run after LINE1 and can overwrite it (no `elif`). In practice L1_ratio and Alu_ratio cannot both be >= 0.8.

### transduction_class (organize_info, L622-706)
- Initialize: `line1_info = None`, `transduction_class = None`
- If repeat_type ends with "LINE1" → `transduction_class = "Solo"`
- If alignment_infos == "---" or is_polyAT is None → early return
- Search LINE1_db for nearby source element:
  - Found → `transduction_class = "Partnered"` (l1_ratio >= 0.01) or `"Orphan"` (l1_ratio < 0.01)
  - Not found → remains "Solo" (for LINE1) or None (for non-LINE1)

### insert_type (annotate_sv_file, L759-771)
Priority order:
1. transduction_class == "Orphan" → Orphan_L1
2. transduction_class == "Partnered" → Partnered_L1
3. transduction_class == "Solo" → Solo_L1
4. repeat_type == "Alu" → Alu
5. repeat_type == "SVA" → SVA
6. ppseudo_info exists → PSD
7. None of the above → "---"

### is_inversion (annotate_sv_file, L773-778)
- repeat_type starts with "Simple" → "Simple"
- repeat_type starts with "Inverted" → "Inverted"
- repeat_type starts with "Other" → "Other"
- Otherwise → "NA"

### LINE1_db search priority (organize_info, L652-706)
Within the tabix fetch window, multiple LINE1_db entries may match. Priority: RepeatMasker > 1000G/HGSVC3 > gnomAD. The `cur_L1_pos` variable selects the nearest source element (closest to alignment position).

## Issues

### 1. L1_SOURCE=None in VCF output (Bug, High)
**Location:** L894 checks `source_info[10]` (transduction_class) to decide whether to emit `L1_SOURCE`. For Solo_L1, transduction_class is "Solo" (passes the check), but `source_info[9]` (line1_info) can be "None" when no source element was found in LINE1_db. Result: `L1_SOURCE=None` in VCF.
**Fix:** Check `source_info[9]` instead:
```python
if source_info[9] not in ("None", "---"):
    annot_info += f";L1_SOURCE={source_info[9]}"
```

### 2. RMSK_INFO semicolons break VCF INFO parsing (Bug, High)
**Location:** `summarize_rmsk()` L357 joins repeat entries with `;`. In VCF, `;` is the INFO field separator, so `RMSK_INFO=entry1;entry2` causes `entry2` to be parsed as a separate INFO tag.
**Fix:** Replace `;` with `|` in VCF output:
```python
if source_info[4] != "---":
    annot_info += f";RMSK_INFO={source_info[4].replace(';', '|')}"
```

### 3. `is_polyAT == None` is always False (Bug, Medium)
**Location:** L628: `if alignment_infos == "---" or is_polyAT == None:`
`is_polyAT` is always a string (from file or default "---"), never Python `None`. This comparison is always False, making it dead code.
**Effect:** When polyA/T is absent (`is_polyAT == "None"`), the early return is skipped. The code proceeds to LINE1_db search with an invalid `source_dir` (defaults to `'-'`). For non-LINE1 insertions this has no practical effect (transduction_class stays None). For LINE1 insertions, it may produce spurious transduction calls.
**Fix:** Change to `is_polyAT in ("None", "---")`:
```python
if alignment_infos == "---" or is_polyAT in ("None", "---"):
```

### 4. IS_INVERSION=NA is redundant for non-LINE1 (Low)
**Location:** VCF output L888 always emits `IS_INVERSION`. For Alu, SVA, PSD, and Unclassified, the value is always "NA" which carries no information.
**Fix (optional):** Only emit for LINE1 types, or omit when value is "NA".

### 5. Python None → string "None" propagation (Low)
**Location:** `organize_info()` uses `str(line1_info)` and `str(transduction_class)` at L709-711. When these are Python `None`, the output is the string `"None"`. Downstream code (`annotate_sv_file`) compares against `"None"` strings, so it works, but is fragile. Using `"---"` consistently as the missing-value sentinel would be cleaner.

### 6. LINE1_db search runs for non-LINE1 when polyA/T is present (Low)
**Location:** organize_info L635-706. If repeat_type is "Alu" or "SVA" but polyA/T is detected, the LINE1_db search still executes. transduction_class stays None (initialized at L624) because it was never set to "Solo". No harm done, but unnecessary work.

## Related files
- `nanomonsv/insert_classify.py`
- `plan/insert_classify_minor.md` (duplicate opens, unused vars, header typos)
