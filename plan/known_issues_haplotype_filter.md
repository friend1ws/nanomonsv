# Known issues: Haplotype Dispersion filter

## Status: WATCH
- Created: 2026-03-16

## Background

Haplotype-aware read counting was implemented in v0.9.0.
Supporting reads are reported per haplotype as `HP1,HP2,HP0` for each breakpoint (BP1, BP2).

A `Haplotype_Dispersion` filter was initially implemented: if `max(HP1,HP2) / (HP1+HP2) < 0.9`, flag as filtered.
The rationale was that somatic SVs should be on one haplotype, so even HP1:HP2 split suggests artifact.

## Problem: filter removes true somatic SVs

Verified with GIAB HG008 somatic benchmark (2026-03-16).

Using threshold `max(HP1,HP2) / (HP1+HP2) < 0.9`:

| Sample | PASS SVs filtered | GIAB somatic matches lost |
|--------|-------------------|--------------------------|
| HG008-T-Revio-GRCh38 | 22 | 4 |
| HG008-T-hiphase-GRCh38 | 26 | 17 (65%) |
| COLO829 | 0 | - |
| H2009 | 0 | - |
| HCC1954 | 0 | - |
| HG008-T-ONT-UL-NE-GRCh38 | 0 | - |
| HG008-T-ONT-UL-NE-CHM13 | 0 | - |

### Why somatic SVs have dispersed haplotype reads

- Somatic SVs are heterozygous in tumor; both WT and mutant alleles present
- Phasing tools assign reads to haplotypes based on nearby germline het variants
- Somatic SV-supporting reads come from the mutant allele, but the mutant allele
  can be on either HP1 or HP2 depending on local phasing
- When tumor purity is high and SV is clonal, most supporting reads are from one HP,
  but some crossover occurs due to:
  - Reads spanning phasing switch errors
  - Reads with ambiguous phasing (assigned to wrong HP)
  - Unphased reads (HP0) diluting the ratio
- hiphase has higher phasing rate → more reads assigned to HP1/HP2 → ratio drops below 0.9 more easily
- ONT-UL-NE has very few phased reads → filter never triggers

### Why cell lines are unaffected

- Cell lines lack matched normal → no germline het calls → no phasing → all reads are HP0
- Filter only activates when HP1+HP2 > 0, so it never fires for unphased data

## Current resolution

- `Haplotype_Dispersion` filter disabled in v0.9.0b2 (post_proc.py calls commented out)
- `--hp_ratio_thres` CLI option commented out in arg_parser.py
- `Haplotype_Dispersion` FILTER header removed from VCF output
- Per-haplotype read counts (HP1,HP2,HP0 for BP1 and BP2) are still reported in output

## Future considerations

- The filter concept is sound for germline artifact removal, but threshold 0.9 is too aggressive for somatic
- A more nuanced approach might work:
  - Only filter when HP1 ≈ HP2 AND HP0 is low (strongly phased data showing true dispersion)
  - Use a lower threshold (e.g., 0.7) combined with minimum phased read count
  - Combine with control panel filtering (germline SVs already removed by panel)
- For now, control panel filtering (max_panel_read_num=2, max_panel_sample_num=2) handles
  germline removal effectively, making the HP filter redundant for that purpose
