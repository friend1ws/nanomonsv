#! /bin/bash

NANOMONSV_PREFIX=$1
REFERENCE_GENOME=$2
SIMPLE_REPEAT_BED=$3

<<_
python3 subscript_postprocess_sbnd/integrate_sbnd.py ${NANOMONSV_PREFIX} ${REFERENCE_GENOME}

python3 subscript_postprocess_sbnd/add_simple_repeat.py ${NANOMONSV_PREFIX}.nanomonsv.proc.result.txt ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp ${SIMPLE_REPEAT_BED}

nanomonsv insert_classify ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt ${REFERENCE_GENOME} --genome_id hg38

python3 subscript_postprocess_sbnd/add_annot_sbnd.py ${NANOMONSV_PREFIX}.nanomonsv.sbnd.proc.result.txt ${SIMPLE_REPEAT_BED} > ${NANOMONSV_PREFIX}.nanomonsv.sbnd.annot.proc.result.txt 
_

Rscript subscript_postprocess_sbnd/plot_sbnd_vis.R ${NANOMONSV_PREFIX} ${NANOMONSV_PREFIX}.nanomonsv.sbnd_vis

rm -rf ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp

