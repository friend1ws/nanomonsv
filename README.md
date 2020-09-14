# nanomonsv

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/friend1ws/nanomonsv.svg?branch=master)](https://travis-ci.org/friend1ws/nanomonsv)

## Introduction

nanomonsv is a software for detecting somatic structural variations from paired (tumor and matched control) cancer genome sequence data. nanomonsv is presented in the following preprint.

Precise characterization of somatic structural variations and mobile element insertions from paired long-read sequencing data with nanomonsv, Shiraishi et al., bioRxiv, 2020, [[link]](https://www.biorxiv.org/content/10.1101/2020.07.22.214262v1)

## Dependency

### For basic use (`parse`, `get` command)

### Binary programs

[htslib](http://www.htslib.org/), [mafft](https://mafft.cbrc.jp/alignment/software/) 

### Python
Pytnon (tested with 3.5, 3.6, 3.7), pysam, numpy, parasail

> [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) (This became optional since version 0.2.0. We have changed the main engine of Smith-Waterman algorithm to parasail.)


### For advanced use (`insert_classify` command)
[bwa](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2), [bedtools](https://bedtools.readthedocs.io/en/latest/), [RepeatMasker](http://www.repeatmasker.org/))

## Preparation

### For basic use (`parse`, `get` command)

#### Install software and add them to the PATH

nanomonsv uses, `tabix`, `bgzip` (which ar part of HTSlib projects) and `mafft` inside the program,
assuming those are installed, and the paths are already added to the running environment.

> ##### For use of SSW Library
> Since version 0.2.0, nanomonsv can be executed without SSW Library. When users want to use SSW Library, create the libssw.so and add the path to the LD_LIBRARY_PATH environment variable. Please refer the **How to use the Python wrapper ssw_lib.py** section in the [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) repository page.


### For advanced use (`insert_classify` command)
`bwa`, `minimap2`, `bedtools` and `RepeatMasker` are required to be installed and these pathese are added to the running environment.


### Input file

nanomonsv accept the BAM file aligned by `minimap2`. 


## Quickstart

1. Install all the prerequisite software and install nanomonsv.
```
pip install nanomonsv (--user)
```

2. Prepare the reference genome for the test data (here, we show the path to Genomic Data Commons reference genome).
```
wget https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834 -O GRCh38.d1.vd1.fa.tar.gz
tar xvf GRCh38.d1.vd1.fa.tar.gz
```

3. Parse the putative structural variation supporting reads of the test data.
```
nanomonsv parse tests/resource/bam/test_tumor.bam output/test_tumor
nanomonsv parse tests/resource/bam/test_ctrl.bam output/test_ctrl
```

4. Get the final result.
```
nanomonsv get output/test_tumor tests/resource/bam/test_tumor.bam GRCh38.d1.vd1.fa --control_prefix output/test_ctrl --control_bam tests/resource/bam/test_ctrl.bam
```

You will see the result file named as `test_tumor.nanomonsv.result.txt`.

## Commands

### parse

This step parses all the supporting reads of putative somatic SVs.

```
nanomonsv parse [-h] [--debug]
                [--split_alignment_check_margin SPLIT_ALIGNMENT_CHECK_MARGIN]
                [--minimum_breakpoint_ambiguity MINIMUM_BREAKPOINT_AMBIGUITY]
                bam_file output_prefix
```
- **bam_file**: Path to input indexed BAM file
- **output_prefix**: Output file prefix

See the help (`nanomonsv parse -h`) for other options.

After successful completion, you will find supporting reads stratified by deletions, insertions, and rearrangements
({output_prefix}.deletion.sorted.bed.gz, {output_prefix}.insertion.sorted.bed.gz, and {output_prefix}.rearrangement.sorted.bedpe.gz)
and their indexes (.tbi files). 


### get

This step gets the SV result from the parsed supporting reads data obtained above.

```
nanomonsv get [-h] [--control_prefix CONTROL_PREFIX]
              [--control_bam CONTROL_BAM]
              [--min_tumor_variant_read_num MIN_TUMOR_VARIANT_READ_NUM]
              [--min_tumor_VAF MIN_TUMOR_VAF]
              [--max_control_variant_read_num MAX_CONTROL_VARIANT_READ_NUM]
              [--max_control_VAF MAX_CONTROL_VAF]
              [--cluster_margin_size CLUSTER_MARGIN_SIZE]
              [--median_mapQ_thres MEDIAN_MAPQ_THRES]
              [--max_overhang_size_thres MAX_OVERHANG_SIZE_THRES]
              [--var_read_min_mapq VAR_READ_MIN_MAPQ] [--debug]
              tumor_prefix tumor_bam reference.fa
 ```
 - **tumor_prefix**: Prefix to the tumor data set in the parse step
 - **tumor_bam**: Path to input indexed BAM file
 - **reference.fa**: Path to reference genome used for the alignment
 
This software can generate a list of SVs without specifying the matched control.
But we have not tested the performance of the approach just using tumor sample, and strongly recommend using the matched control data.
- **control_prefix**: Prefix to the matched control data set in the parse step
- **control_bam**: Path to the matched control BAM file

After successful execution, you will be able to find the result file names as {tumor_prefix}.nanomonsv.result.txt
See the help (`nanomonsv get -h`) for other options. 

When you want to change the engine of Smith-Waterman algorithm to SSW Library, specify `--use_ssw_lib` option,
though we do not generally recomend thi.

#### result

* **Chr_1**: Chromosome for the 1st breakpoint
* **Pos_1**: Coordinate for the 1st breakpoint
* **Dir_1**: Direction of the 1st breakpoint
* **Chr_2**: Chromosome for the 2nd breakpoint
* **Pos_2**: Coordinate for the 2nd breakpoint
* **Dir_2**: Direction of the 2nd breakpoint
* **Inserted_Seq**: Inserted nucleotides within the breakpoints
* **Checked_Read_Num_Tumor**: Total number of reads in the tumor used for the validation alignment step
* **Supporting_Read_Num_Tumor**: Total number of variant reads in the tumor determined in the validation alignment step
* **Checked_Read_Num_Control**: Total number of reads in the normal used for the validation alignment step
* **Supporting_Read_Num_Control**: Total number of variant reads in the matched control determined in the validation alignment step

### insert_classify

This command classifies the long insertions into several mobile element insertions (still in alpha version).

```
nanomonsv insert_classify [-h] [--grc] [--genome_id {hg19,hg38,mm10}]
                          [--debug]
                          sv_list_file output_file reference.fa
```
- **sv_list_file**: SV list file obtained in the get step
- **output_file**: Path to the output file for this command
- **reference.fa**: Path to the reference genome
- **genome_id**: The type of reference genome. Choose from hg19 and hg38 (default is hg38). This is used for selecting LINE1 database.

#### result

* **Insert_Type**: Type of insertion (Solo_L1, Partnered_L1, Orphan_L1, Alu, SVA, PSD)
* **Is_Inversion**: Type of inverted form for Solo LINE1 insertion (Simple, Inverted, Other)
* **L1_Raito**: The match rate with LINE1 sequences for the inserted sequences
* **Alu_Ratio**: The match rate with Alu sequences for the inserted sequences
* **SV_Ratio**: The match rate with SVA sequences for the inserted sequences
* **RMSK_Info**: Summary information of RepeatMasker
* **Alignment_Info**: Alignment information to the human genome
* **Inserted_Pos**: Inserted position (appears only when the inserted sequence is aligned near the inserted positions and implicated to be the tandem duplication or nested LINE1 transduction).
* **Is_PolyA_T**: Extracted poly-A or poly-T sequences
* **Target_Site_Duplication**: Nucleotides of target site duplications
* **L1_Source_Info**: Inferred source site of LINE1 transduction
* **PSD_Gene**: Processed pseudogene name
* **PSD_Overlap_Ratio**: The match rate with the pseudogene
* **PDS_Exon_Num**: The number of pseudogene exons matched with the inserted sequence


### validate

This command, which is part of the procedures of `get` command, 
performs validation of the candidate SVs by alignment of tumor and matched control BAM files.
This may be helpful for the evaluation of SV tools of the short-read platform
when pairs of short-read and long-read sequencing data are available.
This is still in alpha version.

```
nanomonsv validate [-h] [--control_bam CONTROL_BAM]
                   [--var_read_min_mapq VAR_READ_MIN_MAPQ] [--debug]
                   sv_list_file tumor_bam output reference.fa
```
- **sv_list_file**: SV candidate list file (similar format with the result file by `get` command. 
But only from **Chr_1** to **Inserted_Seq** columns are necessary.
- **output_file**: Path to the output file
- **reference.fa**: Path to the reference genome                          
