# Prepare HGSVC3 mobile element insertion database for insert_classify

## Status: TODO
- Created: 2026-03-08

## Overview
Prepare a LINE1 database from the HGSVC3 mobile element insertion (MEI) dataset
for use with `nanomonsv insert_classify`.

Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Mobile_Elements/1.0/

## HGSVC3 MEI data

### Available files (downloaded to `resource/LINE1_db/`)
| File | Description | Records |
|------|-------------|---------|
| `MEI_Callset_GRCh38.ALL.20241211.csv.gz` | Non-reference MEIs called against GRCh38 | 12,642 |
| `MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz` | Non-reference MEIs called against T2T-CHM13 | — |
| `GRCh38_MEIs.ALL.20241211.csv.gz` | Reference MEIs (present in GRCh38 ref, absent in some samples) | 2,468 |
| `T2T-CHM13_MEIs.ALL.20241211.csv.gz` | Reference MEIs for T2T-CHM13 | — |

### MEI Callset format (CSV)
VCF-like columns: ID, CHROM, POS, REF, ALT, QUAL, FILTER, INFO, FORMAT, [sample genotypes...]
Plus: Caller_Count, TE_Designation, L1ME-AID, PALMER, L1ME-AID_INFO, PALMER_INFO, PAVMergedCalls

### TE type breakdown (GRCh38 callset)
| TE_Designation | Count |
|---------------|-------|
| SINE/Alu | 10,270 |
| LINE/L1 | 1,604 |
| Retroposon/SVA | 764 |
| HERVK | 3 |
| snRNA | 1 |

### Key fields for LINE1 db
From PALMER_INFO:
- `SUBFAM` — L1 subfamily (e.g., L1Ta, L1Ambig, L1HS)
- `STRAND` — Insertion orientation (+/-)
- `AF` — Allele frequency
- `SN` — Number of samples with the MEI

From L1ME-AID_INFO:
- `RM_Annotation` — RepeatMasker annotation (e.g., L1HS, AluYa5, SVA_D)
- `Orientation` — Insertion orientation

## Current LINE1 db format

The existing LINE1 db (`nanomonsv/data/LINE1.hg38.bed.gz`) is a tabix-indexed BED file:
```
chr  start  end  info_string  score  strand
```
Where `info_string` = `chr,start,end,strand,name` and sources include:
- RepeatMasker full-length L1 (L1HS, L1PA2-5): 5,944 records
- gnomAD-SV v2.1: 2,546 records
- Total: 8,490 records

The `resource/LINE1_db/LINE1.hg38.bed` is a newer version with gnomAD-SV v3 data (19,254 records).

## Design

### What to extract from HGSVC3
For `insert_classify`, the LINE1 db is used to find source L1 elements near the insertion site
to determine transduction class (Solo vs Partnered vs Orphan).

From HGSVC3 MEI Callset, extract **LINE/L1** records (1,604 in GRCh38) and convert to BED format
compatible with the existing LINE1 db.

### Processing steps
1. Parse `MEI_Callset_GRCh38.ALL.20241211.csv.gz` (and T2T-CHM13 version)
2. Filter for `TE_Designation == "LINE/L1"`
3. Extract: CHROM, POS, POS+1, info_string, AF/SN, STRAND
4. Convert to BED format matching existing LINE1 db
5. Merge with existing sources (RepeatMasker, gnomAD-SV) — avoiding duplicates
6. Sort, bgzip, tabix index

### Open questions
- Should HGSVC3 MEIs be added to the existing LINE1 db, or provided as a separate resource?
- Should Alu and SVA insertions also be included for future use?
- Which reference genome versions to support? (GRCh38, T2T-CHM13, hg19 via liftover?)
- The GRCh38 reference MEIs (`GRCh38_MEIs.ALL.20241211.csv.gz`) represent deletions of
  reference L1 elements — are these useful for insert_classify?

## Related files
- `nanomonsv/insert_classify.py` L529-658 — `organize_info()` uses LINE1_db
- `nanomonsv/data/LINE1.hg38.bed.gz` — Current bundled LINE1 db (8,490 records)
- `resource/LINE1_db/LINE1.hg38.bed` — Newer LINE1 db with gnomAD v3 (19,254 records)
- `resource/LINE1_db/MEI_Callset_GRCh38.ALL.20241211.csv.gz` — Downloaded HGSVC3 data
