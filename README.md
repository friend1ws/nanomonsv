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

See the help (nanomonsv parse -h) for other options.

After successful completion, you will find supporting reads stratified by deletions, insertions, and rearrangements
({output_prefix}.deletion.sorted.bed.gz, {output_prefix}.insertion.sorted.bed.gz, and {output_prefix}.rearrangement.sorted.bedpe.gz)
and their indexes (.tbi files). 



