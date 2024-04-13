#! /usr/bin/env bash

##########
# chain file
if [ ! -f hg38ToHg19.over.chain.gz ]
then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi

if [ ! -f hg38-chm13v2.over.chain.gz ]
then
	wget https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz
fi

if [ ! -f hg19ToHg38.over.chain.gz ]
then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi

if [ ! -f hg19-chm13v2.over.chain.gz ]
then
	wget https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz
fi

##########

##########
# prepare repeat masker file
if [ ! -f rmsk.txt.gz ]
then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
	gzip -f chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
fi

if [ ! -f chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz ]
then
	aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed ./
fi

python proc_rmsk.py rmsk.txt.gz > rmsk.line1.hg38.bed

liftOver rmsk.line1.hg38.bed hg38ToHg19.over.chain.gz rmsk.line1.hg19.bed.tmp rmsk.line1.unmapped

python mod_label.py rmsk.line1.hg19.bed.tmp > rmsk.line1.hg19.bed


python proc_rmsk_chm13.py chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz > rmsk.line1.chm13v2.0.bed 

#########


##########
# 1000 geonme SV file
if [ ! -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz ]
then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
fi

if [ ! -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi ]
then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi
fi

bcftools filter -i 'INFO/SVLEN > 5800 && INFO/SVTYPE == "LINE1"' ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-8 > 1000genomes.line1.hg19.vcf

python proc_1000genomes.py 1000genomes.line1.hg19.vcf > 1000genomes.line1.hg19.bed

liftOver 1000genomes.line1.hg19.bed hg19ToHg38.over.chain.gz 1000genomes.line1.hg38.bed.tmp 1000genomes.line1.hg38.unmapped
liftOver 1000genomes.line1.hg19.bed hg19-chm13v2.over.chain.gz 1000genomes.line1.chm13v2.0.bed.tmp 1000genomes.line1.chm13v2.0.unmapped

python mod_label.py 1000genomes.line1.hg38.bed.tmp > 1000genomes.line1.hg38.bed
python mod_label.py 1000genomes.line1.chm13v2.0.bed.tmp > 1000genomes.line1.chm13v2.0.bed

##########

##########
# gnomAD SV file
if [ ! -f gnomad_v2.1_sv.controls_only.sites.vcf.gz ]
then
    # wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.vcf.gz
    wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.vcf.gz
fi

if [ ! -f gnomad_v2.1_sv.controls_only.sites.vcf.gz.tbi ]
then
    # wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.vcf.gz.tbi
    wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.vcf.gz.tbi
fi

bcftools filter -i 'ALT == "<INS:ME:LINE1>" && INFO/SVLEN >= 5800' gnomad_v2.1_sv.controls_only.sites.vcf.gz | cut -f 1-8 > gnomad.line1.hg19.vcf

python proc_gnomad.py gnomad.line1.hg19.vcf > gnomad.line1.hg19.bed

liftOver gnomad.line1.hg19.bed hg19ToHg38.over.chain.gz gnomad.line1.hg38.bed.tmp gnomad.line1.hg38.unmapped
liftOver gnomad.line1.hg19.bed hg19-chm13v2.over.chain.gz gnomad.line1.chm13v2.0.bed.tmp gnomad.line1.chm13v2.0.unmapped

python mod_label.py gnomad.line1.hg38.bed.tmp > gnomad.line1.hg38.bed
python mod_label.py gnomad.line1.chm13v2.0.bed.tmp > gnomad.line1.chm13v2.0.bed

##########

cat rmsk.line1.hg38.bed 1000genomes.line1.hg38.bed gnomad.line1.hg38.bed | sort -k1,1 -k2,2n -k3,3n > LINE1.hg38.bed

bgzip -c LINE1.hg38.bed > LINE1.hg38.bed.gz

tabix -p bed LINE1.hg38.bed.gz


cat rmsk.line1.hg19.bed 1000genomes.line1.hg19.bed gnomad.line1.hg19.bed | sort -k1,1 -k2,2n -k3,3n > LINE1.hg19.bed

bgzip -c LINE1.hg19.bed > LINE1.hg19.bed.gz

tabix -p bed LINE1.hg19.bed.gz


cat rmsk.line1.chm13v2.0.bed 1000genomes.line1.chm13v2.0.bed gnomad.line1.chm13v2.0.bed | sort -k1,1 -k2,2n -k3,3n > LINE1.chm13v2.0.bed

bgzip -c LINE1.chm13v2.0.bed > LINE1.chm13v2.0.bed.gz

tabix -p bed LINE1.chm13v2.0.bed.gz


mv LINE1.hg38.bed.gz ../nanomonsv/data/
mv LINE1.hg38.bed.gz.tbi ../nanomonsv/data/

mv LINE1.hg19.bed.gz ../nanomonsv/data/
mv LINE1.hg19.bed.gz.tbi ../nanomonsv/data/

mv LINE1.chm13v2.0.bed.gz ../nanomonsv/data/
mv LINE1.chm13v2.0.bed.gz.tbi ../nanomonsv/data/

##########

rm -rf rmsk.line1.hg38.bed 
rm -rf rmsk.line1.hg19.bed
rm -rf rmsk.line1.hg19.bed.tmp 
rm -rf rmsk.line1.unmapped
rm -rf rmsk.line1.chm13v2.0.bed

rm -rf 1000genomes.line1.hg19.vcf
rm -rf 1000genomes.line1.hg19.bed
rm -rf 1000genomes.line1.hg38.bed
rm -rf 1000genomes.line1.hg38.bed.tmp 
rm -rf 1000genomes.line1.hg38.unmapped
rm -rf 1000genomes.line1.chm13v2.0.bed
rm -rf 1000genomes.line1.chm13v2.0.bed.tmp
rm -rf 1000genomes.line1.chm13v2.0.unmapped

rm -rf gnomad.line1.hg19.vcf
rm -rf gnomad.line1.hg38.bed
rm -rf gnomad.line1.hg38.bed.tmp
rm -rf gnomad.line1.hg38.unmapped
rm -rf gnomad.line1.hg19.bed 
rm -rf gnomad.line1.chm13v2.0.bed
rm -rf gnomad.line1.chm13v2.0.bed.tmp
rm -rf gnomad.line1.chm13v2.0.unmapped

rm -rf LINE1.hg38.bed
rm -rf LINE1.hg19.bed
rm -rf LINE1.chm13v2.0.bed

