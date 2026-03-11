# insert_classify improvements (2026-03)

## Status: Done
- Period: 2026-03-08 ~ 2026-03-11

## Summary
Three rounds of improvements to `nanomonsv/insert_classify.py`:
1. VCF input/output support
2. Minor code quality fixes (duplicate opens, unused vars, header typos)
3. Annotation logic audit and sentinel unification

---

## 1. VCF input/output support (2026-03-08 ~ 03-10)

### Overview
Enable `nanomonsv insert_classify` to accept VCF files as input and produce annotated VCF files as output, in addition to the existing TSV-based workflow.

### Design
- **VCF input**: Parse VCF records, extract insertion sequence from SVINSSEQ (INFO) or ALT field. Records with insertion ≥50bp are classified; others pass through as-is.
- **VCF output**: Append INFO fields: INSERT_TYPE, IS_INVERSION, L1_RATIO, ALU_RATIO, SVA_RATIO, RMSK_INFO, TSD, L1_SOURCE, PSD_GENE, PSD_OVERLAP_RATIO, PSD_EXON_NUM.
- **Format detection**: `.vcf` or `.vcf.gz` → VCF mode, otherwise → TSV mode.
- **Classification threshold**: Changed from 100bp to 50bp (both modes).

### Changes
- Added `extract_insert_seq_from_vcf_record()`, `make_fasta_file_from_vcf()`, `annotate_vcf_file()`, `classify_record()`
- Modified `insert_classify_main()` in `run.py` for format detection

---

## 2. Minor code quality fixes (2026-03-09 ~ 03-10)

### Fixes applied
- **Duplicate file opens**: Removed unused `pysam.AlignmentFile` and `open()` calls in `summarize_bwa_alignment` and `summarize_bwa_alignment2`
- **Unused variables**: Removed `is_secondary`/`is_supplementary` in `sam2bed_split`
- **Header typos**: Fixed `L1_Raito`→`L1_Ratio`, `SV_Ratio`→`SVA_Ratio`, `PDS_Exon_Num`→`PSD_Exon_Num`

---

## 3. Annotation logic audit and refactoring (2026-03-10 ~ 03-11)

### Data flow
```
inserted sequence
    ├── RepeatMasker  → summarize_rmsk()           → repeat_type, L1/Alu/SVA_ratio, rmsk_info
    ├── BWA alignment → summarize_bwa_alignment2()  → alignment_infos, inserted_pos
    └── check_tsd_polyAT()                          → is_polyAT, tsd
            │
            ▼
    organize_info()  ← LINE1_db
        → source_info[0..10]:
            [0] repeat_type        "Simple_LINE1", "Inverted_LINE1", "Other_LINE1", "Alu", "SVA", "---"
            [1] l1_ratio           float as str
            [2] alu_ratio          float as str
            [3] sva_ratio          float as str
            [4] rmsk_info          semicolon-separated repeat entries
            [5] alignment_infos    semicolon-separated alignment entries
            [6] inserted_pos       "chr,start,end" or "---"
            [7] is_polyAT          "polyA", "polyT", or "---"
            [8] tsd                sequence or "---"
            [9] line1_info         LINE1 source ID or "---"
           [10] transduction_class "Solo", "Partnered", "Orphan", or "---"
            │
            ▼
    annotate_sv_file() / annotate_vcf_file()
        → insert_type:   Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD, "---"
        → is_inversion:  Simple, Inverted, Other, NA
        + ppseudo_info:  PSD_Gene, PSD_Overlap_Ratio, PSD_Exon_Num
```

Missing values are represented as `"---"` throughout.

### Bugs fixed
1. **L1_SOURCE=None in VCF** — Was checking `source_info[10]` (transduction_class) instead of `source_info[9]` (line1_info). Fixed to check `source_info[9] != "---"`.
2. **RMSK_INFO semicolons break VCF** — `;` in RMSK_INFO conflicted with VCF INFO separator. Fixed by replacing `;` with `|` in VCF output.
3. **`is_polyAT == None` always False** — `is_polyAT` is a string from file, never Python `None`. Removed dead condition from early return.
4. **None → "None" propagation** — Unified all missing-value sentinels to `"---"` across `check_tsd_polyAT`, `proc_rmsk_info`, `organize_info`, `annotate_sv_file`.

### Refactoring
- Extracted `_load_annotation_data()` shared by `annotate_sv_file` and `annotate_vcf_file`
- Simplified LINE1_db priority logic: replaced 3 booleans (`is_rmsk`, `is_1000g`, `is_gnomad`) with single `is_ref_source`. ref_source (rmsk) takes priority over non-ref (HGSVC3, 1000genomes, gnomAD) which are treated equally.
- Removed unused `overlap_ratio` variable

### Deferred (Low priority)
- IS_INVERSION=NA redundant for non-LINE1 types
- LINE1_db search runs unnecessarily for non-LINE1 with polyA/T

---

## Related files
- `nanomonsv/insert_classify.py`
- `nanomonsv/run.py` — `insert_classify_main()`
- `nanomonsv/arg_parser.py` — Argument parser for insert_classify
