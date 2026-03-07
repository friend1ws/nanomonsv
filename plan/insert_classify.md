# insert_classify: VCF input/output support

## Status: TODO
- Created: 2026-03-08

## Overview
Enable `nanomonsv insert_classify` to accept VCF files as input and produce
annotated VCF files as output, in addition to the current TSV-based workflow.

## Current state

### Input format (TSV)
`insert_classify` currently expects the nanomonsv result.txt (TSV), where:
- Columns: Chr_1, Pos_1, Dir_1, Chr_2, Pos_2, Dir_2, Inserted_Seq, ...
- `make_fasta_file()` extracts `Inserted_Seq` (column index 6) for insertions ≥100bp
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
- Parse VCF to extract insertion records (SVTYPE=INS or SVTYPE=BND with SVINSSEQ)
- Extract SVINSSEQ from INFO field as the insertion sequence (equivalent to Inserted_Seq)
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
- Consider adding `--output_format {tsv,vcf}` option, or auto-detect from output extension

### Open questions
- Should the VCF output be a new file, or should annotations be added in-place to the input VCF?
- For BND records with SVINSSEQ: should these also be classified?
  (Currently, insert_classify only processes insertions with Dir_1='+' and Dir_2='-')
- How to handle VCF records without SVINSSEQ (e.g., symbolic `<INS>` without sequence)?

## Related files
- `nanomonsv/insert_classify.py` — All classification functions
- `nanomonsv/run.py` L509-586 — `insert_classify_main()`
- `nanomonsv/arg_parser.py` L211-238 — Argument parser for insert_classify
- `nanomonsv/vcf_convert.py` — Existing VCF conversion (no insert_classify integration)
