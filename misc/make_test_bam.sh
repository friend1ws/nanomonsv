#! /usr/bin/env bash

if [ ! -d test_bam_tmp ]
then
    mkdir test_bam_tmp
fi

TUMOR_BAM=s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/GDC.GRCh38.d1.vd1/minimap-2.17_genebay-201903_guppy-3.4.5/COLO829/COLO829.bam

samtools view -bh ${TUMOR_BAM} chr10:87940338-87940738 > test_bam_tmp/t1.bam
samtools view -bh ${TUMOR_BAM} chr10:87952384-87952784 > test_bam_tmp/t2.bam
samtools view -bh ${TUMOR_BAM} chr11:62010282-62010683 > test_bam_tmp/t3.bam
samtools view -bh ${TUMOR_BAM} chr12:129287032-129287435 > test_bam_tmp/t4.bam
samtools view -bh ${TUMOR_BAM} chr15:84141774-84142174 > test_bam_tmp/t5.bam
samtools view -bh ${TUMOR_BAM} chr7:151049371-151049771 > test_bam_tmp/t6.bam
samtools view -bh ${TUMOR_BAM} chr18:68712023-68712423 > test_bam_tmp/t7.bam
samtools view -bh ${TUMOR_BAM} chr18:68715389-68715789 > test_bam_tmp/t8.bam
samtools view -bh ${TUMOR_BAM} chrX:78265187-78265587 > test_bam_tmp/t9.bam

samtools cat test_bam_tmp/t1.bam test_bam_tmp/t2.bam test_bam_tmp/t3.bam test_bam_tmp/t4.bam test_bam_tmp/t5.bam test_bam_tmp/t6.bam test_bam_tmp/t7.bam test_bam_tmp/t8.bam test_bam_tmp/t9.bam > test_bam_tmp/t.bam
samtools sort test_bam_tmp/t.bam > test_tumor.bam
samtools index test_tumor.bam
rm -rf test_bam_tmp/*

CONTROL_BAM=s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/GDC.GRCh38.d1.vd1/minimap-2.17_genebay-201903_guppy-3.4.5/COLO829BL/COLO829BL.bam

samtools view -bh ${CONTROL_BAM} chr10:87940338-87940738 > test_bam_tmp/t1.bam
samtools view -bh ${CONTROL_BAM} chr10:87952384-87952784 > test_bam_tmp/t2.bam
samtools view -bh ${CONTROL_BAM} chr11:62010282-62010683 > test_bam_tmp/t3.bam
samtools view -bh ${CONTROL_BAM} chr12:129287032-129287435 > test_bam_tmp/t4.bam
samtools view -bh ${CONTROL_BAM} chr15:84141774-84142174 > test_bam_tmp/t5.bam
samtools view -bh ${CONTROL_BAM} chr7:151049371-151049771 > test_bam_tmp/t6.bam
samtools view -bh ${CONTROL_BAM} chr18:68712023-68712423 > test_bam_tmp/t7.bam
samtools view -bh ${CONTROL_BAM} chr18:68715389-68715789 > test_bam_tmp/t8.bam
samtools view -bh ${CONTROL_BAM} chrX:78265187-78265587 > test_bam_tmp/t9.bam

samtools cat test_bam_tmp/t1.bam test_bam_tmp/t2.bam test_bam_tmp/t3.bam test_bam_tmp/t4.bam test_bam_tmp/t5.bam test_bam_tmp/t6.bam test_bam_tmp/t7.bam test_bam_tmp/t8.bam test_bam_tmp/t9.bam > test_bam_tmp/t.bam
samtools sort test_bam_tmp/t.bam > test_ctrl.bam
samtools index test_ctrl.bam
rm -rf test_bam_tmp
