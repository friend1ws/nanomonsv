# nanomonsv

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CI](https://github.com/friend1ws/nanomonsv/actions/workflows/python-test.yml/badge.svg)

## Introduction

nanomonsv is a software for detecting somatic structural variations from paired (tumor and matched control) cancer genome sequence data. nanomonsv is presented in the following paper. **When you use nanomonsv or any resource of this repository, please kindly cite this paper**.

Precise characterization of somatic complex structural variations from tumor/control paired long-read sequencing data with nanomonsv, Shiraishi et al., Nucleic Acids Research, 2023, [[link]](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad526/7201946).

Key features:
- **Single-nucleotide breakpoint resolution** using consensus sequences from long-read alignments.
- **LINE1 insertion classification**: Distinguishes Solo L1, Partnered L1 (transduction), and Orphan L1 (orphan transduction), and identifies source L1 elements.
- **Two detection modules**: Canonical SV module for standard SVs with high precision and recall, and [Single breakend SV module](https://github.com/friend1ws/nanomonsv/wiki/Single-breakend-SV) for complex SVs involving highly-repetitive sequences (centromeres, LINE1, viruses) that can only be identified by long-reads.

## Installation

```
pip install nanomonsv
```

You can also install via conda (bioconda channel). Occasionally the conda releases lag behind PyPI.
```
conda create -n nanomonsv -c conda-forge -c bioconda nanomonsv
```

### Dependencies

| Tool | Required for | Notes |
|------|-------------|-------|
| [htslib](http://www.htslib.org/) (tabix, bgzip) | parse, get | Must be in PATH |
| [mafft](https://mafft.cbrc.jp/alignment/software/) | get | Consensus generation |
| [racon](https://github.com/isovic/racon) | get (`--use_racon`) | Recommended |
| [bwa](https://github.com/lh3/bwa) | insert_classify | |
| [minimap2](https://github.com/lh3/minimap2) | insert_classify | |
| [bedtools](https://bedtools.readthedocs.io/en/latest/) | insert_classify | |
| [RepeatMasker](http://www.repeatmasker.org/) | insert_classify | |

Python >=3.9, pysam, numpy, parasail

### Input requirements

- BAM or CRAM file aligned by minimap2
- For CRAM files, specify `--reference_fasta`


## Quick Start

1. Prepare the reference genome (here, GDC reference genome).
```
wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz
tar xvf GRCh38.d1.vd1.fa.tar.gz
```

2. Parse putative SV supporting reads.
```
nanomonsv parse tests/resource/bam/test_tumor.bam output/test_tumor
nanomonsv parse tests/resource/bam/test_ctrl.bam output/test_ctrl
```

3. Get the final result.
```
nanomonsv get output/test_tumor tests/resource/bam/test_tumor.bam GRCh38.d1.vd1.fa \
    --control_prefix output/test_ctrl --control_bam tests/resource/bam/test_ctrl.bam
```

You will find the result file `test_tumor.nanomonsv.result.txt`.

## Usage

### parse

Parses all supporting reads of putative somatic SVs.

```
nanomonsv parse [-h] [--reference_fasta reference.fa] [--debug]
                [--split_alignment_check_margin SPLIT_ALIGNMENT_CHECK_MARGIN]
                [--minimum_breakpoint_ambiguity MINIMUM_BREAKPOINT_AMBIGUITY]
                alignment_file output_prefix
```
- **alignment_file**: Path to input indexed BAM or CRAM file
- **output_prefix**: Output file prefix
- **--reference_fasta**: Path to reference genome (recommended for CRAM files)

After successful completion, you will find:
`{output_prefix}.deletion.sorted.bed.gz`, `{output_prefix}.insertion.sorted.bed.gz`, `{output_prefix}.rearrangement.sorted.bedpe.gz`, `{output_prefix}.bp_info.sorted.bed.gz` and their indexes (.tbi files).


### get

Gets the SV result from parsed supporting reads.

```
nanomonsv get [-h] [--control_prefix CONTROL_PREFIX]
              [--control_bam CONTROL_BAM]
              [--control_panel_prefix CONTROL_PANEL_PREFIX]
              [--simple_repeat_bed SIMPLE_REPEAT_BED]
              [--min_tumor_variant_read_num MIN_TUMOR_VARIANT_READ_NUM]
              [--min_tumor_VAF MIN_TUMOR_VAF]
              [--max_control_variant_read_num MAX_CONTROL_VARIANT_READ_NUM]
              [--max_control_VAF MAX_CONTROL_VAF]
              [--cluster_margin_size CLUSTER_MARGIN_SIZE]
              [--median_mapQ_thres MEDIAN_MAPQ_THRES]
              [--max_overhang_size_thres MAX_OVERHANG_SIZE_THRES]
              [--var_read_min_mapq VAR_READ_MIN_MAPQ]
              [--qv10] [--qv15] [--qv20] [--qv25] [--use_racon]
              [--single_bnd] [--processes PROCESSES]
              [--sort_option SORT_OPTION] [--max_memory_minimap2] [--debug]
              tumor_prefix tumor_bam reference.fa
```
- **tumor_prefix**: Prefix to the tumor data set in the parse step
- **tumor_bam**: Path to input indexed BAM file
- **reference.fa**: Path to reference genome used for the alignment

#### Recommended options

| Option | Recommendation | Description |
|--------|---------------|-------------|
| `--control_prefix` / `--control_bam` | **Strongly recommended** | Matched control for somatic filtering. We strongly recommend using matched control data whenever possible. |
| `--control_panel_prefix` | Recommended | Non-matched control panel (see [Control Panel](#control-panel)) |
| `--simple_repeat_bed` | **Strongly recommended** | Filter indels in simple repeats. BED files provided in [resource/simple_repeats](https://github.com/friend1ws/nanomonsv/tree/master/resource/simple_repeats) |
| `--use_racon` | Recommended | Better breakpoint resolution and polished inserted sequences |
| `--single_bnd` | Optional | Enable single breakend SV detection (use with `--use_racon`). See [wiki](https://github.com/friend1ws/nanomonsv/wiki/Single-breakend-SV) |
| `--processes N` | Optional | Multi-processing mode |

#### Quality presets

| Preset | Recommended for |
|--------|----------------|
| `--qv10` | ONT data with median Q10 (e.g., Guppy before v5) |
| `--qv15` | ONT data with median Q15 (e.g., Guppy v5/v6) |
| `--qv20` | ONT data with median Q20+ (e.g., Dorado SUP, Q20+ chemistry) |
| `--qv25` | PacBio HiFi data |

### merge_control

Merges non-matched control panel supporting reads obtained by `parse`.

```
nanomonsv merge_control [-h] prefix_list_file output_prefix
```

- **prefix_list_file**: List of output_prefix generated at the `parse` stage
- **output_prefix**: Prefix to the merged control supporting reads


### insert_classify

Classifies long insertions into mobile element insertions (LINE1, Alu, SVA, processed pseudogene).

```
nanomonsv insert_classify [-h] [--debug] sv_list_file output_file reference.fa gencode.gtf.gz LINE1_db
```
- **sv_list_file**: SV list file obtained in the get step
- **output_file**: Path to the output file
- **reference.fa**: Path to the reference genome
- **gencode.gtf.gz**: Path to gene annotation GTF file. We recommend [Gencode basic annotation](https://www.gencodegenes.org/human/) (e.g., [gencode.v49.basic.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz))
- **LINE1_db**: Path to LINE1 database. Use the files in [resource/LINE1_db](https://github.com/friend1ws/nanomonsv/tree/master/resource/LINE1_db)


### validate

Validates candidate SVs by alignment of tumor and matched control BAM files.
This may be helpful for evaluating SV tools on short-read platforms
when pairs of short-read and long-read sequencing data are available.

```
nanomonsv validate [-h] [--control_bam CONTROL_BAM]
                   [--var_read_min_mapq VAR_READ_MIN_MAPQ] [--debug]
                   sv_list_file tumor_bam output reference.fa
```
- **sv_list_file**: SV candidate list file (only **Chr_1** to **Inserted_Seq** columns are necessary)
- **output_file**: Path to the output file
- **reference.fa**: Path to the reference genome


## Output Format

### Canonical SV result ({tumor_prefix}.nanomonsv.result.txt)

| Column | Description |
|--------|------------|
| Chr_1 | Chromosome for the 1st breakpoint |
| Pos_1 | Coordinate for the 1st breakpoint |
| Dir_1 | Direction of the 1st breakpoint |
| Chr_2 | Chromosome for the 2nd breakpoint |
| Pos_2 | Coordinate for the 2nd breakpoint |
| Dir_2 | Direction of the 2nd breakpoint |
| Inserted_Seq | Inserted nucleotides within the breakpoints (`---` if none) |
| SV_ID | Identifier of SVs |
| Checked_Read_Num_Tumor | Total reads in the tumor used for validation alignment |
| Supporting_Read_Num_Tumor | Variant reads in the tumor from validation alignment |
| Checked_Read_Num_Control | Total reads in the matched control used for validation alignment |
| Supporting_Read_Num_Control | Variant reads in the matched control from validation alignment |
| Is_Filter | Filter status (`PASS` or filter reason such as `Simple_repeat`) |

A VCF format file ({tumor_prefix}.nanomonsv.result.vcf) is also generated.
See the [wiki page](https://github.com/friend1ws/nanomonsv/wiki/How-to-understand-nanomonsv-result-filtering) for details on filtering.

### Single breakend result ({tumor_prefix}.nanomonsv.sbnd.result.txt)

Generated when `--single_bnd` is specified.

| Column | Description |
|--------|------------|
| Chr_1 | Chromosome of the breakpoint |
| Pos_1 | Coordinate of the breakpoint |
| Dir_1 | Direction of the breakpoint |
| Contig | Assembled contig sequence at the breakpoint |
| SV_ID | Identifier of the single breakend |
| Checked_Read_Num_Tumor | Total reads in the tumor used for validation alignment |
| Supporting_Read_Num_Tumor | Variant reads in the tumor from validation alignment |
| Checked_Read_Num_Control | Total reads in the matched control used for validation alignment |
| Supporting_Read_Num_Control | Variant reads in the matched control from validation alignment |
| Is_Filter | Filter status (`PASS`, `Simple_repeat`, `Canonical_SV_overlap`, or combinations) |

A VCF format file ({tumor_prefix}.nanomonsv.sbnd.result.vcf) is also generated,
using VCF single breakend notation (e.g., `N.` or `.N` in ALT field with `SVTYPE=BND`).

### insert_classify result

| Column | Description |
|--------|------------|
| Insert_Type | Type of insertion (Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD) |
| Is_Inversion | Inverted form for Solo LINE1 (Simple, Inverted, Other) |
| L1_Ratio | Match rate with LINE1 sequences |
| Alu_Ratio | Match rate with Alu sequences |
| SVA_Ratio | Match rate with SVA sequences |
| RMSK_Info | Summary information of RepeatMasker |
| Alignment_Info | Alignment information to the human genome |
| Inserted_Pos | Inserted position (for tandem duplication or nested LINE1 transduction) |
| Is_PolyA_T | Extracted poly-A or poly-T sequences |
| Target_Site_Duplication | Nucleotides of target site duplications |
| L1_Source_Info | Inferred source site of LINE1 transduction |
| PSD_Gene | Processed pseudogene name |
| PSD_Overlap_Ratio | Match rate with the pseudogene |
| PSD_Exon_Num | Number of pseudogene exons matched with the inserted sequence |


## Control Panel

We recommend using a control panel for filtering common SVs and noise.
Use `merge_control` to create one from your own sequencing data.

For users without sufficient control data, we provide a pre-built control panel
created from 30 Nanopore sequencing samples from the [Human Pangenome Reference Consortium](https://humanpangenome.org/),
available at [zenodo](https://zenodo.org/records/11470934).

This control panel was built by aligning 30-40 Nanopore (basecalled by Guppy ver. 4 and 6) and PacBio HiFi sequencing data
to the GRCh38/CHM13 reference genome with minimap2 version 2.24.
We recommend using a control panel as close as possible in platform and basecall quality.
When unsure, use the Nanopore control panel from Guppy version 4 (noisier panels tend to be more versatile).

**When you use these control panels and publish, do not forget to credit [HPRC](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#citations)!**


## Example Data

The Oxford Nanopore Sequencing data used in the paper is available through the public sequence repository (BioProject ID: PRJDB10898):
- COLO829: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248304[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248305[accn])
- H2009: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248308[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248309[accn])
- HCC1954: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248306[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248307[accn])

Results of nanomonsv for the above data are available [here](https://github.com/friend1ws/nanomonsv/tree/master/misc/example).
Please kindly cite the [NAR paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad526/7201946) when you use these data.

See the tutorial [wiki page](https://github.com/friend1ws/nanomonsv/wiki/Tutorial) for an example workflow on analyzing the COLO829 sample.


## Citation

Shiraishi et al., Precise characterization of somatic complex structural variations from tumor/control paired long-read sequencing data with nanomonsv, Nucleic Acids Research, 2023, [[link]](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad526/7201946).
