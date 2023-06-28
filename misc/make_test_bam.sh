#! /usr/bin/env bash

set -eux

if [ ! -d test_bam_tmp ]
then
    mkdir test_bam_tmp
fi

TUMOR_BAM=s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/GDC.GRCh38.d1.vd1/minimap-2.17_genebay-201903_guppy-3.4.5/COLO829/COLO829.bam

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr1:224458701-224459101 > test_bam_tmp/t01.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr1:224612218-224612618 > test_bam_tmp/t02.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr10:7017350-7017750    > test_bam_tmp/t03.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr10:7090710-7091110    > test_bam_tmp/t04.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr10:58717264-58717861  > test_bam_tmp/t05.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr12:72272912-72273497  > test_bam_tmp/t06.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr15:23440220-23440620  > test_bam_tmp/t07.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr15:84141773-84142173  > test_bam_tmp/t08.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr19:17285802-17286202  > test_bam_tmp/t09.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr19:17286631-17287031  > test_bam_tmp/t11.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr20:15019776-15020176  > test_bam_tmp/t12.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr20:15032996-15033396  > test_bam_tmp/t13.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr20:38645624-38646024  > test_bam_tmp/t14.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr3:25358911-25359311   > test_bam_tmp/t15.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr3:25359368-25359768   > test_bam_tmp/t16.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr3:26390229-26390629   > test_bam_tmp/t17.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr6:26193611-26194011   > test_bam_tmp/t18.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${TUMOR_BAM} chr7:151049370-151049770 > test_bam_tmp/t19.bam

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools cat  test_bam_tmp/t01.bam test_bam_tmp/t02.bam test_bam_tmp/t03.bam test_bam_tmp/t04.bam test_bam_tmp/t05.bam \
  test_bam_tmp/t06.bam test_bam_tmp/t07.bam test_bam_tmp/t08.bam test_bam_tmp/t09.bam test_bam_tmp/t11.bam \
  test_bam_tmp/t12.bam test_bam_tmp/t13.bam test_bam_tmp/t14.bam test_bam_tmp/t15.bam test_bam_tmp/t16.bam \
  test_bam_tmp/t17.bam test_bam_tmp/t18.bam test_bam_tmp/t19.bam > test_bam_tmp/t.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools sort test_bam_tmp/t.bam > test_tumor.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools index test_tumor.bam

rm -rf test_bam_tmp/*

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -hC test_tumor.bam -T ~/resources/database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa > test_tumor.cram
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools index test_tumor.cram

CONTROL_BAM=s3://eva-bucket-tokyo/kataoka-lab/long_read_sequencing/cell-line/GDC.GRCh38.d1.vd1/minimap-2.17_genebay-201903_guppy-3.4.5/COLO829BL/COLO829BL.bam

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr1:224458701-224459101 > test_bam_tmp/t01.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr1:224612218-224612618 > test_bam_tmp/t02.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr10:7017350-7017750    > test_bam_tmp/t03.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr10:7090710-7091110    > test_bam_tmp/t04.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr10:58717264-58717861  > test_bam_tmp/t05.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr12:72272912-72273497  > test_bam_tmp/t06.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr15:23440220-23440620  > test_bam_tmp/t07.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr15:84141773-84142173  > test_bam_tmp/t08.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr19:17285802-17286202  > test_bam_tmp/t09.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr19:17286631-17287031  > test_bam_tmp/t11.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr20:15019776-15020176  > test_bam_tmp/t12.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr20:15032996-15033396  > test_bam_tmp/t13.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr20:38645624-38646024  > test_bam_tmp/t14.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr3:25358911-25359311   > test_bam_tmp/t15.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr3:25359368-25359768   > test_bam_tmp/t16.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr3:26390229-26390629   > test_bam_tmp/t17.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr6:26193611-26194011   > test_bam_tmp/t18.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -bh ${CONTROL_BAM} chr7:151049370-151049770 > test_bam_tmp/t19.bam

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools cat  test_bam_tmp/t01.bam test_bam_tmp/t02.bam test_bam_tmp/t03.bam test_bam_tmp/t04.bam test_bam_tmp/t05.bam \
  test_bam_tmp/t06.bam test_bam_tmp/t07.bam test_bam_tmp/t08.bam test_bam_tmp/t09.bam test_bam_tmp/t11.bam \
  test_bam_tmp/t12.bam test_bam_tmp/t13.bam test_bam_tmp/t14.bam test_bam_tmp/t15.bam test_bam_tmp/t16.bam \
  test_bam_tmp/t17.bam test_bam_tmp/t18.bam test_bam_tmp/t19.bam > test_bam_tmp/t.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools sort test_bam_tmp/t.bam > test_ctrl.bam
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools index test_ctrl.bam

rm -rf test_bam_tmp

/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools view -hC test_ctrl.bam -T ~/resources/database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa > test_ctrl.cram
/home/aiokada/conda/x64/envs/nanomonsv_devel/bin/samtools index test_ctrl.cram
