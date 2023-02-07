#! /bin/bash

NANOMONSV_PREFIX=$1

python3 subscript_postprocess_sbnd/integrate_sbnd.py ${NANOMONSV_PREFIX} ~/environment/data/reference/GRCh38.d1.vd1.fa

python3 subscript_postprocess_sbnd/add_simple_repeat.py ${NANOMONSV_PREFIX}.nanomonsv.proc.result.txt ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp simpleRepeat.bed.gz   

nanomonsv subscript_postprocess_sbnd/insert_classify ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt ~/environment/data/reference/GRCh38.d1.vd1.fa --genome_id hg38

python3 subscript_postprocess_sbnd/add_annot_sbnd.py ${NANOMONSV_PREFIX}.nanomonsv.sbnd.proc.result.txt simpleRepeat.bed.gz > ${NANOMONSV_PREFIX}.nanomonsv.sbnd.annot.proc.result.txt 

Rscript subscript_postprocess_sbnd/plot_sbnd_vis.R ${NANOMONSV_PREFIX} ${NANOMONSV_PREFIX}.nanomonsv.sbnd_vis

rm -rf ${NANOMONSV_PREFIX}.nanomonsv.annot.proc.result.txt.tmp

