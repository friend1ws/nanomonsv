# Roadmap

Release plan from v0.8.2 (current) to v0.9.0.

## v0.8.2 — Bug fix + sbnd improvements + sbnd VCF (Released)

Low risk. No change to existing canonical SV output format.

- [x] [vcf_fix](done/vcf_fix.md) — Fix BND ALT field bug (tref2→tref1), sort order, missing FILTER header
- [x] [sbnd_simple_repeat](done/sbnd_simple_repeat.md) — Simple repeat filtering for sbnd (only when `--simple_repeat_bed` is set)
- [x] [sbnd_dedup_canonical](done/sbnd_dedup_canonical.md) — Remove sbnd overlapping canonical SVs (coordinate comparison only, no RepeatMasker/BWA dependency)
- [x] [sbnd_vcf](done/sbnd_vcf.md) — VCF output for single breakends (VCF 4.3 single breakend format with contig sequence)

## v0.9.0 — insert_classify extension + misc cleanup + default options

- [x] [insert_classify_improve](done/insert_classify_improve2603.md) — VCF input/output, minor fixes, annotation logic audit & sentinel unification
- [x] [hgsvc3_mei_db](done/hgsvc3_mei_db.md) — Prepare HGSVC3 MEI database (LINE1 from 1,604 records + Alu/SVA)
- [x] [misc_cleanup](done/misc_cleanup.md) — Remove misc/ scripts whose functionality has been integrated into `nanomonsv get`
- [x] [readme_improve](readme_improve.md) — Key features section, GTF description, quality presets
- [ ] [default_options](default_options.md) — Make racon and `--single_bnd` default on in `nanomonsv get`

## v0.9.0 bugfix — GIAB workflow test issues

- [ ] [giab_bugfix](giab_bugfix.md) — Fix crash on long sbnd contigs (csv field limit), cap consensus length at 100KB, suppress minimap2 core dumps

## v0.9.1 — Haplotype support

- [ ] [haplotag](haplotag.md) — Per-haplotype supporting read counts (HP1/HP2/unphased), haplotype dispersion filter

## v0.9.2 — Micro homology

- [ ] [micro_homology](micro_homology.md) — Detect micro homology at breakpoints, add HOMLEN/HOMSEQ/CIPOS to VCF

## Dependency graph

```
v0.8.2  vcf_fix                 ✓
        sbnd_simple_repeat      ✓
        sbnd_dedup_canonical    ✓
        sbnd_vcf                ✓
            │
            ▼
v0.9.0  insert_classify_improve ✓ + hgsvc3_mei_db ✓ + misc_cleanup ✓
        + readme_improve ✓ + default_options
            │
            ▼
v0.9.1  haplotag
            │
            ▼
v0.9.2  micro_homology
```
