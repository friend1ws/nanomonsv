# Roadmap

Release plan from v0.8.1 (current) to v0.9.0.

## v0.8.2 — Bug fix + sbnd improvements

Low risk. No change to existing canonical SV output format.

- [x] [vcf_fix](done/vcf_fix.md) — Fix BND ALT field bug (tref2→tref1), sort order, missing FILTER header
- [ ] [sbnd_simple_repeat](sbnd_simple_repeat.md) — Simple repeat filtering for sbnd (only when `--simple_repeat_bed` is set)
- [ ] [sbnd_dedup_canonical](sbnd_dedup_canonical.md) — Remove sbnd overlapping canonical SVs (coordinate comparison only, no RepeatMasker/BWA dependency)

Rationale: The VCF bug is critical and should be fixed immediately. The two sbnd
changes are self-contained within `post_proc.py` and do not affect canonical SV
output. Easy to regression-test with existing tests.

## v0.8.3 — VCF output improvements

Changes to VCF output fields and breakpoint reporting.

- [ ] [sbnd_vcf](sbnd_vcf.md) — VCF output for single breakends (VCF 4.3 single breakend format with contig sequence)
- [ ] [micro_homology](micro_homology.md) — Detect micro homology at breakpoints, add HOMLEN/HOMSEQ/CIPOS to VCF

Rationale: Both are VCF output changes, efficient to test together. micro_homology
affects breakpoint position interpretation, so it should be separated from the
bug fix release. sbnd_vcf benefits from the sbnd quality improvements in v0.8.2.

## v0.8.4 — Haplotype support

Changes to intermediate file format and output columns.

- [ ] [haplotag](haplotag.md) — Per-haplotype supporting read counts (HP1/HP2/unphased), haplotype dispersion filter

Rationale: Modifies intermediate file format in `count_sread_by_alignment.py` and
adds output columns to result.txt. Needs careful testing to ensure backward
compatibility when BAM has no HP tags. Best released independently.

## v0.9.0 — insert_classify extension

New functionality for mobile element insertion classification.

- [ ] [insert_classify](insert_classify.md) — VCF input/output support for insert_classify
- [ ] [hgsvc3_mei_db](hgsvc3_mei_db.md) — Prepare HGSVC3 MEI database (LINE1 from 1,604 records + Alu/SVA)

Rationale: VCF I/O for insert_classify is a significant feature addition, warranting
a minor version bump. The HGSVC3 database update is tightly coupled with
insert_classify and should be released together.

## Dependency graph

```
v0.8.2  vcf_fix
        sbnd_simple_repeat
        sbnd_dedup_canonical
            │
            ▼
v0.8.3  sbnd_vcf  (after sbnd quality improves)
        micro_homology
            │
            ▼
v0.8.4  haplotag  (independent, needs thorough testing)
            │
            ▼
v0.9.0  insert_classify + hgsvc3_mei_db
```
