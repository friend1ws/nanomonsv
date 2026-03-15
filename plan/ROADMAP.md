# Roadmap

Release plan from v0.8.2 (current) to v0.9.0.

## v0.8.2 — Bug fix + sbnd improvements + sbnd VCF (Released)

Low risk. No change to existing canonical SV output format.

- [x] [vcf_fix](done/vcf_fix.md) — Fix BND ALT field bug (tref2→tref1), sort order, missing FILTER header
- [x] [sbnd_simple_repeat](done/sbnd_simple_repeat.md) — Simple repeat filtering for sbnd (only when `--simple_repeat_bed` is set)
- [x] [sbnd_dedup_canonical](done/sbnd_dedup_canonical.md) — Remove sbnd overlapping canonical SVs (coordinate comparison only, no RepeatMasker/BWA dependency)
- [x] [sbnd_vcf](done/sbnd_vcf.md) — VCF output for single breakends (VCF 4.3 single breakend format with contig sequence)

## v0.9.0 — insert_classify extension + misc cleanup + default options (Released)

- [x] [insert_classify_improve](done/insert_classify_improve2603.md) — VCF input/output, minor fixes, annotation logic audit & sentinel unification
- [x] [hgsvc3_mei_db](done/hgsvc3_mei_db.md) — Prepare HGSVC3 MEI database (LINE1 from 1,604 records + Alu/SVA)
- [x] [misc_cleanup](done/misc_cleanup.md) — Remove misc/ scripts whose functionality has been integrated into `nanomonsv get`
- [~] [readme_improve](readme_improve.md) — Key features section, GTF description, quality presets (mostly done; Zenodo URL, panel threshold docs, get option list pending)
- [x] [default_options](done/default_options.md) — Make racon and `--single_bnd` default on in `nanomonsv get`
- [x] [haplotag](done/haplotag.md) — Per-haplotype supporting read counts (HP1/HP2/unphased)
- [x] [giab_bugfix](done/giab_bugfix.md) — Fix crash on long sbnd contigs (csv field limit), minimap2 `-f 0.05` for satellite repeats, insert_classify negative start fix

## v0.9.0 post-release — Control panel tuning + known issues

- [x] Default panel thresholds: `max_panel_read_num=2`, `max_panel_sample_num=2` (tuned for 1000-sample ONT panels)
- [x] Haplotype_Dispersion filter disabled (removes true somatic SVs; see [known_issues_haplotype_filter](known_issues_haplotype_filter.md))
- [x] 1kg-ont-vienna control panels: singleton removal for size reduction, Zenodo upload
- [ ] [known_issues_generate_consensus](known_issues_generate_consensus.md) — sbnd consensus length cap, RLIMIT_CORE suppression (mitigated, not fully resolved)
- [ ] [known_issues_haplotype_filter](known_issues_haplotype_filter.md) — HP filter too aggressive for somatic SVs (disabled, future improvement needed)

## Future

- [ ] [insert_classify_2604](insert_classify_2604.md) — Documentation, literature review, label/naming review, TXT/VCF compatibility
- [ ] [micro_homology](micro_homology.md) — Detect micro homology at breakpoints, add HOMLEN/HOMSEQ/CIPOS to VCF
- [x] [fix_zero_length_exon_in_sam2bed_split](done/fix_zero_length_exon_in_sam2bed_split.md) — Fix zero-length exon bug in sam2bed_split

## Dependency graph

```
v0.8.2  vcf_fix                 ✓
        sbnd_simple_repeat      ✓
        sbnd_dedup_canonical    ✓
        sbnd_vcf                ✓
            │
            ▼
v0.9.0  insert_classify_improve ✓ + hgsvc3_mei_db ✓ + misc_cleanup ✓
        + readme_improve ✓ + default_options ✓
        + haplotag ✓ + giab_bugfix ✓
            │
            ▼
Future  insert_classify_2604
        micro_homology
        fix_zero_length_exon_in_sam2bed_split
```
