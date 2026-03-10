# README.md improvement

## Status: TODO
- Created: 2026-03-08

## Overview
Restructure README.md to improve readability without losing essential information.

## Current problems

### 1. Flat structure
All commands (parse, get, merge_control, insert_classify, validate) are at the same heading level.
The primary workflow (parse → get) is not visually distinguished from advanced/optional commands.

### 2. get section is overloaded
The `get` section mixes:
- Command syntax and arguments
- Tips/recommendations (control, control panel, simple repeat, racon, preset)
- Result format description

Tips take up more space than the actual command documentation.

### 3. Redundancy between Dependency and Preparation
Both sections list the same tools (htslib, mafft, racon, bwa, minimap2, etc.)
with slightly different framing. Can be merged.

### 4. Scattered version history
Phrases like "From v0.3.0", "From v0.4.0", "From v0.5.0", "From v0.7.0", "From v0.8.0"
appear throughout. These are useful as changelogs but clutter the main documentation.

### 5. HTML tags in markdown
`<ins>` tags are used for underlining subsection headers. Non-standard and renders
inconsistently across markdown viewers.

### 6. Missing sbnd documentation
sbnd result format and VCF output are not documented (to be added in v0.8.2 release).

## Proposed structure

```
# nanomonsv
  [badges]
  ## Introduction
  ## Installation
    - pip / conda
    - Dependencies (merged from Dependency + Preparation)
    - Input file requirements
  ## Quick Start
    (existing Quickstart, trimmed)
  ## Usage
    ### parse
    ### get
      - Basic usage
      - Options
        (control, control panel, simple repeat, racon, preset — as a table or compact list)
    ### merge_control
    ### insert_classify
    ### validate
  ## Output Format
    ### Canonical SV result (.nanomonsv.result.txt)
    ### Canonical SV VCF (.nanomonsv.result.vcf)
    ### Single breakend result (.nanomonsv.sbnd.result.txt)  [NEW in v0.8.2]
    ### Single breakend VCF (.nanomonsv.sbnd.result.vcf)  [NEW in v0.8.2]
  ## Control Panel
    (moved from inside get section — important enough for its own section)
  ## Example Data
    (existing Realistic example data)
  ## Citation
```

## Specific changes

### Merge Dependency + Preparation → Installation
Current:
```
## Dependency
### For basic use
### Binary programs
htslib, mafft, racon
### Python
Python, pysam, numpy, parasail
### For advanced use
bwa, minimap2, bedtools, RepeatMasker

## Preparation
### For basic use
#### Install software and add them to the PATH
(repeats htslib, mafft info)
###### For use of racon
(repeats racon info)
### For advanced use
(repeats bwa, minimap2, bedtools, RepeatMasker)
### Input file
### Control panel
```

Proposed:
```
## Installation
pip install nanomonsv
conda install -c bioconda nanomonsv

### Dependencies
| Tool | Required for | Notes |
|------|-------------|-------|
| htslib (tabix, bgzip) | parse, get | Must be in PATH |
| mafft | get | Consensus generation |
| racon | get (--use_racon) | Recommended |
| bwa | insert_classify | |
| minimap2 | insert_classify | |
| bedtools | insert_classify | |
| RepeatMasker | insert_classify | |

### Input requirements
- BAM/CRAM aligned by minimap2
- For CRAM: use --reference_fasta
```

### Compact get options
Replace verbose Tips subsections with a recommendation table:

| Option | Recommendation | Description |
|--------|---------------|-------------|
| --control_prefix/--control_bam | Strongly recommended | Matched control for filtering |
| --control_panel_prefix | Recommended | Non-matched control panel |
| --simple_repeat_bed | Strongly recommended | Filter indels in simple repeats |
| --use_racon | Recommended | Better breakpoint resolution |
| --single_bnd | Optional | Detect single breakend SVs |
| --qv10/15/20/25 | Match your data | Quality preset |

Keep brief explanations below the table for items that need more context
(e.g., control panel download link, simple repeat BED file location).

### Remove or relocate version history notes
- "From v0.3.0, we support racon" → just document racon as an option
- "From v0.4.0, VCF format" → just document VCF in Output Format
- "From v0.5.0, control panel" → just document control panel
- "From v0.7.0, preset/CRAM" → just document preset and CRAM support
- "From v0.8.0, simple_repeat_bed" → just document the option
- "For the older versions (before v0.4.0)" → remove (link to wiki if needed)

Optionally add a brief "What's New" section at the top or link to CHANGELOG.

### Separate Output Format section
Move result column descriptions out of `get` section into dedicated Output Format section.
Add sbnd result and VCF descriptions.

### Replace `<ins>` tags
Use bold (`**text**`) or `###` subheadings instead.

## Files to modify
- `README.md`

## Notes
- Do at release time (v0.8.2) together with sbnd documentation
- Keep wiki links intact — wiki pages have additional details
- Keep citation, example data, and control panel download links unchanged
