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

### Preparation

### Install softwares and add them to the PATH

nanomonsv uses, tabix, bgzip (which ar part of HTSlib projects) and mafft inside the program,
assuming those are installed and the pathes are already added to the running environment.
Also, 
