# nanomonsv

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/friend1ws/nanomonsv.svg?branch=master)](https://travis-ci.com/friend1ws/nanomonsv.svg?branch=master&status=started)

## Introduction

nanomonsv is a software for detecting somatic structural variations from paired (tumor and matched control) cancer genome sequence data. nanomonsv is presented in the following preprint. **When you use nanomonsv or any resource of this repository, please kindly site this preprint**.

Precise characterization of somatic complex structural variations from paired long-read sequencing data with nanomonsv, Shiraishi et al., bioRxiv, 2020, [[link]](https://www.biorxiv.org/content/10.1101/2020.07.22.214262v3).

The current version of nanomonsv includes two detection modules, Canonical SV module, and [Single breakend SV module](https://github.com/friend1ws/nanomonsv/wiki/Single-breakend-SV). Canonical SV module can identify somatic SVs that can be captured by short-read technologies with higher precision and recall than existing methods. Furthermore, Single breakend SV module enables the detection of complex SVs that can only be identified by long-reads, such as SVs involving highly-repetitive centromeric sequences, and LINE1- and virus-mediated rearrangements. 

Please see the [wiki page]([Single breakend SV module](https://github.com/friend1ws/nanomonsv/wiki/Single-breakend-SV)) for Single breakend SV module.

## Dependency

### For basic use (`parse`, `get` command)

### Binary programs

[htslib](http://www.htslib.org/), [mafft](https://mafft.cbrc.jp/alignment/software/), [racon](https://github.com/isovic/racon)(optional from ver. 0.3.0. However, we recommend to use this option. Add --use_racon option when you perfrom get command.)

### Python
Pytnon (tested with >=3.6), pysam, numpy, parasail

> [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) (This became optional since version 0.2.0. We have changed the main engine of Smith-Waterman algorithm to parasail.)


### For advanced use (`insert_classify` command)
[bwa](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2), [bedtools](https://bedtools.readthedocs.io/en/latest/), [RepeatMasker](http://www.repeatmasker.org/)

## Preparation

### For basic use (`parse`, `get` command)

#### Install software and add them to the PATH

nanomonsv uses, `tabix`, `bgzip` (which ar part of HTSlib projects) and `mafft` inside the program,
assuming those are installed, and the paths are already added to the running environment.

> ##### For use of SSW Library
> Since version 0.2.0, nanomonsv can be executed without SSW Library. When users want to use SSW Library, create the libssw.so and add the path to the LD_LIBRARY_PATH environment variable. Please refer the **How to use the Python wrapper ssw_lib.py** section in the [SSW Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) repository page.

###### For use of racon
Since version 0.3.0, we support racon for the step where generating consensus sequence and get single-base resolution breakpoints. racon may become the default instead of mafft in the future.


### For advanced use (`insert_classify` command)
`bwa`, `minimap2`, `bedtools` and `RepeatMasker` are required to be installed and these pathese are added to the running environment.


### Input file

nanomonsv accept the BAM file aligned by `minimap2`. 


### Control panel
Starting with version 0.5.0, the use of the control panel is supported. 
In this software, supporting reads for SVs are collected for multiple samples other than the target sample, 
and such reads are removed as common noise (or those derived from common SVs) in the `get` stage. 
This strategy is expected to exclude many false positives as well as improve computational cost.

We have prepared the command to create control panels from the user's own sequencing data.
In addition, for users who do not have sufficient sequencing data that can serve as a control panel (or just do not have time for processing), 
we prepared a control panel that has been created using the 30 Nanopore sequencing data from the [Human Pangenome Reference Consortium](https://humanpangenome.org/), which is available at [zenodo](https://zenodo.org/record/7017953).

This control panel is made by aligning 30 Nanopore sequencing data to the GRCh38 reference genome (obtained from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects)) with minimap2 version 2.24. 
**When you use these control panels and publish, do not forget to credit to [HPRC](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#citations)!**


## Quickstart

1. Install all the prerequisite software and install nanomonsv.
```
pip install nanomonsv (--user)
```

You can also install nanomonsv via conda (bioconda channel).
```
conda create -n nanomonsv -c conda-forge -c bioconda nanomonsv
```
Occasionally the conda releases lag behind the source code and PyPI releases.

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

You will see the result file named `test_tumor.nanomonsv.result.txt`.

## Realistic example sequencing data

The Oxford Nanopore Sequencing data used in the bioRxiv paper is available through the public sequence repository service (BioProject ID: PRJDB10898):
- COLO829: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248304[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248305[accn])
- H2009: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248308[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248309[accn])
- HCC1954: [tumor](https://www.ncbi.nlm.nih.gov/sra/DRX248306[accn]), [control](https://www.ncbi.nlm.nih.gov/sra/DRX248307[accn])

The results of nanomonsv for the above data are available [here](https://github.com/friend1ws/nanomonsv/tree/master/misc/example).
When you perform nanomonsv to the above data and have experienced errors, please report to us.
Also, please kindly cite the [bioRxiv paper](https://www.biorxiv.org/content/10.1101/2020.07.22.214262v1) when you use these data.

See tutorial [wiki page](https://github.com/friend1ws/nanomonsv/wiki/Tutorial) for an example workflow on analyzing COLO829 sample.

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

After successful completion, you will find supporting reads stratified by deletions, insertions, and rearrangements: 
({output_prefix}.deletion.sorted.bed.gz, {output_prefix}.insertion.sorted.bed.gz, {output_prefix}.rearrangement.sorted.bedpe.gz, and {output_prefix}.bp_info.sorted.bed.gz and their indexes (.tbi files). 


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
              [--var_read_min_mapq VAR_READ_MIN_MAPQ] [--use_ssw_lib] [--use_racon]
              [--single_bnd] [--threads THREADS] [--processes PROCESSES] 
              [--sort_option SORT_OPTION] [--max_memory_minimap2] [--debug]
              tumor_prefix tumor_bam reference.fa
 ```
 - **tumor_prefix**: Prefix to the tumor data set in the parse step
 - **tumor_bam**: Path to input indexed BAM file
 - **reference.fa**: Path to reference genome used for the alignment
 
This software can generate a list of SVs without specifying the matched control.
But we have not tested the performance of the approach just using tumor sample, and strongly recommend using the matched control data.
Even when only tumor sample is available, we still recommend using dummy control sample (collected from other person's tissue).
- **control_prefix**: Prefix to the matched control data set in the parse step
- **control_bam**: Path to the matched control BAM file

When you use the control panel (recommended!), use the following argument.
- **control_panel_prefix**: Prefix of non-matched control panel data processed in merge_control step.

After successful execution, you will be able to find the result file names as {tumor_prefix}.nanomonsv.result.txt.
See the help (`nanomonsv get -h`) for other options. 

When you want to change the engine of Smith-Waterman algorithm to SSW Library, specify `--use_ssw_lib` option,
though we do not generally recommend this.

Also, we basically recommend using `--use_racon` option. This will slightly improve the identification of single-base resolution breakpoint, 
and polishing of inserted sequences. 

For detection of single breakend SVs, please use `--single_bnd` option as well as `--use_racon`. 
Please wee [wiki page](https://github.com/friend1ws/nanomonsv/wiki/Single-breakend-SV).

Also, we have prepared the script (misc/post_fileter.py) for filtering the result.
Please see the [wiki page](https://github.com/friend1ws/nanomonsv/wiki/How-to-filter-nanomonsv-result).
For output files of the version 0.4.0 and later, some filtering has already been performed (see the [wiki page](https://github.com/friend1ws/nanomonsv/wiki/How-to-understand-nanomonsv-result-filtering)). 
However, we strongly recommed to perform additional processing; removing indels within simple repeat regions (see the [wiki page](https://github.com/friend1ws/nanomonsv/wiki/An-example-on-removing-indels-within-simple-repeat)).

From the version 0.4.0, we will also provide the VCF format result files.

#### result

* **Chr_1**: Chromosome for the 1st breakpoint
* **Pos_1**: Coordinate for the 1st breakpoint
* **Dir_1**: Direction of the 1st breakpoint
* **Chr_2**: Chromosome for the 2nd breakpoint
* **Pos_2**: Coordinate for the 2nd breakpoint
* **Dir_2**: Direction of the 2nd breakpoint
* **Inserted_Seq**: Inserted nucleotides within the breakpoints
* **SV_ID**: Identifier of SVs (originally comming from cluster of SV supporting reads)
* **Checked_Read_Num_Tumor**: Total number of reads in the tumor used for the validation alignment step
* **Supporting_Read_Num_Tumor**: Total number of variant reads in the tumor determined in the validation alignment step
* **Checked_Read_Num_Control**: Total number of reads in the normal used for the validation alignment step
* **Supporting_Read_Num_Control**: Total number of variant reads in the matched control determined in the validation alignment step
* **Is_Filter**: Filter status. PASS if this SV has passed all the filters


### merge_control

This command merges non-matched control panel supporting reads obtained by performing `parse` command.

```
nanomonsv merge_control [-h] prefix_list_file output_prefix
```

- **prefix_list_file**: The list of output_prefix generated at the above `parse` stage. 
- **output_prefix**: Prefix to the merged control supporting reads. 


### insert_classify

This command classifies the long insertions into several mobile element insertions (still in alpha version).
This does not yet support VCF format, but we will do so in the near future.

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
* **SVA_Ratio**: The match rate with SVA sequences for the inserted sequences
* **RMSK_Info**: Summary information of RepeatMasker
* **Alignment_Info**: Alignment information to the human genome
* **Inserted_Pos**: Inserted position (appears only when the inserted sequence is aligned near the other insertion and implicated to be the tandem duplication or nested LINE1 transduction).
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
