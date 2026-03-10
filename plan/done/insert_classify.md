# insert_classify: VCF input/output support

## Status: DONE
- Created: 2026-03-08
- Completed: 2026-03-10

## Overview
Enable `nanomonsv insert_classify` to accept VCF files as input and produce
annotated VCF files as output, in addition to the current TSV-based workflow.

## Current state

### Input format (TSV)
`insert_classify` currently expects the nanomonsv result.txt (TSV), where:
- Columns: Chr_1, Pos_1, Dir_1, Chr_2, Pos_2, Dir_2, Inserted_Seq, ...
- `make_fasta_file()` extracts `Inserted_Seq` (column index 6) for insertions ≥100bp → **change to ≥50bp**
- The key for each record is `chr1,pos1,dir1,chr2,pos2,dir2,len(inserted_seq)`

### Output format (TSV)
`annotate_sv_file()` appends 14 columns to the original TSV:
- Insert_Type (Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD, ---)
- Is_Inversion (Simple, Inverted, Other, NA)
- L1_Ratio, Alu_Ratio, SV_Ratio
- RMSK_Info, Alignment_Info, Inserted_Pos
- Is_PolyA_T, Target_Site_Duplication, L1_Source_Info
- PSD_Gene, PSD_Overlap_Ratio, PDS_Exon_Num

### Processing pipeline (`insert_classify_main()`)
1. `make_fasta_file()` — Extract insertion sequences to FASTA
2. Processed pseudogene detection (minimap2 + bedtools + exon GTF)
3. RepeatMasker — Classify repeat content of insertions
4. `check_tsd_polyAT()` — Detect target site duplications and poly-A/T tails
5. BWA alignment — Align insertions to reference genome
6. `organize_info()` — Combine all annotations, determine L1 transduction class
7. `annotate_sv_file()` — Merge annotations back into original TSV

### Relationship with VCF conversion
- `nanomonsv get` outputs both result.txt (TSV) and result.vcf
- `genomesv2vcf_convert()` runs on result.txt **before** insert_classify
- insert_classify results are **not reflected in VCF** — only in the annotated TSV

## Design

### Approach
Support VCF as both input and output for `insert_classify`.

#### VCF input
- Parse all VCF records; extract insertion sequence from each:
  - SVINSSEQ in INFO field (preferred)
  - Otherwise, ALT field sequence (remove REF prefix/suffix)
- Records with insertion sequence ≥50bp are classified; others pass through as-is
- No SVTYPE filtering — any record with a long enough insertion is a target
- Convert VCF records to the internal key format used by the pipeline
- The rest of the processing pipeline (RepeatMasker, BWA, etc.) remains unchanged

#### VCF output
- After classification, write results back as VCF with additional INFO fields:
  - `INSERT_TYPE` — Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD
  - `IS_INVERSION` — Simple, Inverted, Other
  - `RMSK_INFO` — RepeatMasker classification summary
  - `TSD` — Target site duplication sequence
  - `L1_SOURCE` — Source L1 element info for transductions
  - `PSD_GENE` — Processed pseudogene source gene
- Add corresponding `##INFO` header lines

### Input detection
Detect input format by file extension or content:
- `.vcf` or `.vcf.gz` → VCF mode
- Otherwise → TSV mode (current behavior, backward compatible)

### Changes required

#### 1. `nanomonsv/insert_classify.py`
- Add `make_fasta_file_from_vcf()` — VCF equivalent of `make_fasta_file()`
- Add `annotate_vcf_file()` — VCF equivalent of `annotate_sv_file()`
- Modify `insert_classify_main()` (in `run.py`) to detect format and call appropriate functions

#### 2. `nanomonsv/arg_parser.py`
- The `sv_list_file` argument already accepts a path; no change needed for VCF input
- Auto-detect output format from output file extension (`.vcf` → VCF, otherwise → TSV)

### Decisions
- **Classification threshold**: 50bp (changed from 100bp). Applies to both TSV and VCF modes.
- **Classification target**: Any record with insertion sequence ≥50bp, regardless of SVTYPE.
- **Sequence extraction (VCF)**: SVINSSEQ (INFO) preferred; if absent, extract from ALT field (`ALT[len(REF):]` for simple INS where `len(ALT) > len(REF)` and both are base sequences). BND ALT notation (e.g., `A[chr2:pos2[TCGATCG`) is not supported.
- **Internal key format**: VCF records use the same 7-element key as TSV (`chr,pos,+,chr,pos,-,len`) for pipeline compatibility. `CHROM_POS_ID` is used only for matching results back to VCF records at output.
- **VCF output**: New file (not in-place modification of input VCF)
- **`<INS>` without sequence**: Skip with warning message to stderr (e.g., "Skipping record {ID}: symbolic `<INS>` without sequence is not supported")
- **Non-target records in VCF output**: Pass through as-is (no classification INFO added). "Non-target" = no insertion sequence or <50bp.

## Related files
- `nanomonsv/insert_classify.py` — All classification functions
- `nanomonsv/run.py` L509-586 — `insert_classify_main()`
- `nanomonsv/arg_parser.py` L211-238 — Argument parser for insert_classify
- `nanomonsv/vcf_convert.py` — Existing VCF conversion (no insert_classify integration)
