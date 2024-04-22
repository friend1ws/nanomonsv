#!/bin/bash

mkdir ../workspace
wget -O ../workspace/human_GRCh37_simpleRepeat.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz   
zcat ../workspace/human_GRCh37_simpleRepeat.txt.gz  | cut -f 2-4 | sort -k1,1 -k2,2n -k3,3n > ../workspace/human_GRCh37_simpleRepeat.bed   
bgzip -c ../workspace/human_GRCh37_simpleRepeat.bed > ../human_GRCh37_simpleRepeat.bed.gz
tabix -p bed ../human_GRCh37_simpleRepeat.bed.gz

wget -O ../workspace/human_GRCh38_simpleRepeat.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz   
zcat ../workspace/human_GRCh38_simpleRepeat.txt.gz | cut -f 2-4 | sort -k1,1 -k2,2n -k3,3n > ../workspace/human_GRCh38_simpleRepeat.bed   
bgzip -c ../workspace/human_GRCh38_simpleRepeat.bed > ../human_GRCh38_simpleRepeat.bed.gz
tabix -p bed ../human_GRCh38_simpleRepeat.bed.gz

wget -O ../workspace/human_chm13_repeatMasker.bed.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
grep "Simple_repeat" ../workspace/human_chm13_repeatMasker.bed.gz | cut -f 1-3 | sort -k 1,1 -k 2,2n -k 3,3n > ../workspace/human_chm13v2.0_simpleRepeat.bed
bgzip -c ../workspace/human_chm13v2.0_simpleRepeat.bed > ../human_chm13v2.0_simpleRepeat.bed.gz
tabix -p bed ../human_chm13v2.0_simpleRepeat.bed.gz
