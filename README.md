# nanomonsv

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/friend1ws/nanomonsv.svg?branch=master)](https://travis-ci.org/friend1ws/nanomonsv)

## Introduction

nanomonsv is a software for detecting somatic structural variass from paired (tumor and matched control) cancer genome sequence data. 

## Dependency

### Binary programs
[htslib](http://www.htslib.org/), [mafft](https://mafft.cbrc.jp/alignment/software/), [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

### Python
Pytnon (>= 3.6), pysam, numpy, scipy, statistics, swalign

## Preparation

### Install softwares and add them to the PATH

nanomonsv uses, tabix, bgzip (which ar part of HTSlib projects) and mafft inside the program,
assuming those are installed and the pathes are already added to the running environment.
Also, for the preparation of SSW Library, 
create the libssw.so and add the path to the LD_LIBRARY_PATH environment variable.

### Input file

nanomonsv the input file aligned by `minimap2`. 


## Commands

### parse

This step parse all the supporting reads of putative somatic SVs.

```
nanomonsv parse [-h] [--debug]
                [--split_alignment_check_margin SPLIT_ALIGNMENT_CHECK_MARGIN]
                [--minimum_breakpoint_ambiguity MINIMUM_BREAKPOINT_AMBIGUITY]
                bam_file output_prefix
```
- **bam_file**: Path to input indexed bam file
- **output_prefix**: Output file prefix

See the help (`nanomonsv parse -h`) for other options.

After successful completion, you will find supporting reads stratified by deletions, insertions, and rearrangements
({output_prefix}.deletion.sorted.bed.gz, {output_prefix}.insertion.sorted.bed.gz, and {output_prefix}.rearrangement.sorted.bedpe.gz)
and their indexes (.tbi files). 


### get

This step get the SV result from the parsed supporting reads data obtained above.

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
 - **tumor_bam**: Path to input indexed bam file
 - **reference.fa**: Path to reference genome used for the alignment
 
This software can generate the list of SVs without specifying the matched control.
But we have not tested the performance of the approach just using tumor sample, and recommend to use the matched control data.
- **control_prefix**: Prefix to the matched control data set set in the parse step
- **control_bam**: Path to the matched control bam file

After successful execution, you will be able to find the result file names as {tumor_prefix}.nanomonsv.result.txt
See the help (`nanomonsv get -h`) for other options. 


### insert_classify

This step classify the long insertions into several mobile element insertions (still in alpha version).

```
nanomonsv insert_classify [-h] [--grc] [--genome_id {hg19,hg38,mm10}]
                          [--debug]
                          sv_list_file output_file reference.fa
```
- **sv_list_file**: SV list file obtained in the get step
- **output_file**: Path to the output file for this command
- **reference.fa**: Path to the reference genome


