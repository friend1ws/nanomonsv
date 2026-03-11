# insert_classify: documentation & annotation improvement

## Status: TODO
- Created: 2026-03-11

## Overview
Improve understanding, documentation, and annotation quality of `insert_classify`.

## Tasks

### 1. Developer documentation: internal logic reference
Write a detailed internal reference for the insert_classify pipeline.
Destination: `docs/insert_classify_internals.md` (or similar)

Should cover:
- Full pipeline flow with intermediate files (.tmp.fasta, .tmp.rmsk.txt, .tmp.bwa.sam, .tmp.tsd.polyAT.txt, .tmp.alignment.txt, .tmp.org.txt, .tmp.ppseudo.txt)
- Each function's role, input/output format, and key decisions
- source_info index mapping and what each field means concretely
- Alignment_Info format: `qstart,qend,strand,mapQ,mismatch,ins,del,chr,rstart,rend,p/s` — what each field is
- Inserted_Pos: derived from Alignment_Info, indicates whether insertion maps near the original SV breakpoint (±5000bp)
- LINE1_db search logic: direction selection (source_dir), search window, priority (is_ref_source)
- transduction_class determination: Solo vs Partnered vs Orphan
- Classification priority in annotate_sv_file/annotate_vcf_file

### 2. User documentation: README / wiki
Add user-facing explanation of insert_classify results.
Destination: README.md (insert_classify section) and/or wiki page

Should cover:
- VCF input/output support (currently undocumented)
- Insert_Type meanings:
  - Solo_L1: Insertion contains mostly LINE1 sequence (L1_Ratio ≥ 0.8), no transduction detected
  - Partnered_L1: LINE1 + flanking sequence from source locus (transduction), source identified in LINE1_db
  - Orphan_L1: Flanking transduced sequence without LINE1 body (L1_Ratio < 0.01), source identified
  - Alu: Alu element insertion (Alu_Ratio ≥ 0.8)
  - SVA: SVA element insertion (SVA_Ratio ≥ 0.8)
  - PSD: Processed pseudogene insertion (multi-exon match to known transcript)
- Is_Inversion: whether LINE1 has 5' inversion (common structural feature of L1 insertions)
- TSD: target site duplication — hallmark of retrotransposition
- Is_PolyA_T: poly-A/T tail detection (hallmark of LINE1/Alu insertion mechanism)
- L1_Source_Info: inferred source LINE1 element for transductions
- VCF INFO fields (currently 3 fields missing compared to TSV: Alignment_Info, Inserted_Pos, Is_PolyA_T)

### 3. Literature review: what annotations should be reported
Review existing tools and papers to understand best practices for MEI annotation.

#### Tools to review
- **TraFiC** (Tubio et al., Science 2017) — LINE1 transduction detection in cancer
  - What output fields does it report?
  - Source element identification method
  - Transduction classification (solo, partnered/TD+, orphan/TD0)
- **MELT** (Gardner et al., Genome Research 2017) — MEI detection from short reads
  - Output VCF INFO fields for MEIs
- **xTea** (Chu et al., Nature Communications 2021) — Transposable element analysis
  - Classification categories and output format
- **GRAFFITE** — Long-read MEI genotyping
- **Palmer** (Ewing, Mobile DNA 2019) — Reference LINE1 source annotation

#### Key questions
- What additional annotations are commonly reported?
  - Source element genomic coordinates (chr, start, end, strand)
  - Source element name/family (L1HS, L1PA2, etc.)
  - Transduced sequence length
  - 5' truncation point
  - EN cleavage site motif
  - Internal inversion breakpoint position
  - Allele frequency of source element (from population databases)
- Which of these can nanomonsv already compute (or derive from existing intermediate data)?
- Which would require new analysis steps?

### 4. Label/naming review
Review and potentially rename output fields for clarity and consistency.

Current concerns:
- `Is_Inversion`: "Simple/Inverted/Other/NA" — "Simple" means non-inverted, which is confusing
  - Consider: "None/5prime_inversion/Complex/NA" or similar
- `Insert_Type` vs `transduction_class`: internal naming inconsistency
  - `transduction_class` = "Solo/Partnered/Orphan" (organize_info)
  - `Insert_Type` = "Solo_L1/Partnered_L1/Orphan_L1/Alu/SVA/PSD" (annotate)
  - The "_L1" suffix is added at annotation time — is this clear enough?
- `L1_Source_Info`: field name suggests it's just L1, but format includes coordinates + db entry name
  - Consider splitting into separate fields: source_chr, source_pos, source_strand, source_name
- `RMSK_Info`: compact format hard to parse — should VCF use a more structured representation?
- `Alignment_Info`: not in VCF output — should it be? Or is it too detailed for VCF?
- `Inserted_Pos`: name is ambiguous — "inserted position" could mean where the insertion is, not where it maps to. Consider "Mapped_Pos" or "Source_Alignment_Pos"

## Priority
1. Task 3 (literature review) — informs decisions for all other tasks
2. Task 4 (label review) — should be settled before writing documentation
3. Task 1 (developer docs)
4. Task 2 (user docs + README update)
