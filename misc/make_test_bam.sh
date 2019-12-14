#! /usr/bin/env bash

if [ ! -d test_bam_tmp ]
then
    mkdir test_bam_tmp
fi

samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 10:89700199-89700399 > test_bam_tmp/t1.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 10:89712241-89712441 > test_bam_tmp/t2.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 7:150746557-150746757 > test_bam_tmp/t3.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 15:84810625-84810825 > test_bam_tmp/t4.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 18:66379362-66379562 > test_bam_tmp/t5.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829.bam 18:66382728-66382928 > test_bam_tmp/t6.bam
samtools cat test_bam_tmp/t1.bam test_bam_tmp/t2.bam test_bam_tmp/t3.bam test_bam_tmp/t4.bam test_bam_tmp/t5.bam test_bam_tmp/t6.bam > test_bam_tmp/t.bam
samtools sort test_bam_tmp/t.bam > test_tumor.bam
samtools index test_tumor.bam
rm -rf test_bam_tmp/*

samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 10:89700199-89700399 > test_bam_tmp/t1.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 10:89712241-89712441 > test_bam_tmp/t2.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 7:150746557-150746757 > test_bam_tmp/t3.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 15:84810625-84810825 > test_bam_tmp/t4.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 18:66379362-66379562 > test_bam_tmp/t5.bam
samtools view -bh s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/minimap2/COLO829BL.bam 18:66382728-66382928 > test_bam_tmp/t6.bam
samtools cat test_bam_tmp/t1.bam test_bam_tmp/t2.bam test_bam_tmp/t3.bam test_bam_tmp/t4.bam test_bam_tmp/t5.bam test_bam_tmp/t6.bam > test_bam_tmp/t.bam
samtools sort test_bam_tmp/t.bam > test_ctrl.bam
samtools index test_ctrl.bam
rm -rf test_bam_tmp
