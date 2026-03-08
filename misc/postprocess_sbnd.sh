#! /bin/bash
set -e

NANOMONSV_PREFIX=$1
REFERENCE_GENOME=$2

python3 subscript_sbnd/annotate_contig.py ${NANOMONSV_PREFIX} ${REFERENCE_GENOME}

Rscript subscript_sbnd/plot_contig.R ${NANOMONSV_PREFIX} ${NANOMONSV_PREFIX}.nanomonsv.sbnd_vis
