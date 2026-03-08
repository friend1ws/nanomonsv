#! /usr/bin/env bash
set -euo pipefail

##########
# Download chain files
for chain in \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
    https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
    https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz
do
    fname=$(basename "$chain")
    [ -f "$fname" ] || wget "$chain"
done

##########
# RepeatMasker
[ -f rmsk.txt.gz ] || wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

if [ ! -f chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz ]; then
    [ -f chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed ] || \
        aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed ./
    sort -k1,1 -k2,2n -k3,3n chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed | bgzip -c > chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz
    tabix -p bed chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz
fi

python subscript/proc_rmsk.py rmsk.txt.gz > rmsk.line1.hg38.bed
liftOver rmsk.line1.hg38.bed hg38ToHg19.over.chain.gz rmsk.line1.hg19.bed.tmp rmsk.line1.unmapped
python subscript/mod_label.py rmsk.line1.hg19.bed.tmp > rmsk.line1.hg19.bed

python subscript/proc_rmsk_chm13.py chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz > rmsk.line1.chm13v2.0.bed

##########
# 1000 Genomes SV
[ -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz ] || \
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
[ -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi ] || \
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi

bcftools filter -i 'INFO/SVLEN > 5800 && INFO/SVTYPE == "LINE1"' ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-8 > 1000genomes.line1.hg19.vcf
python subscript/proc_1000genomes.py 1000genomes.line1.hg19.vcf > 1000genomes.line1.hg19.bed

liftOver 1000genomes.line1.hg19.bed hg19ToHg38.over.chain.gz 1000genomes.line1.hg38.bed.tmp 1000genomes.line1.hg38.unmapped
liftOver 1000genomes.line1.hg19.bed hg19-chm13v2.over.chain.gz 1000genomes.line1.chm13v2.0.bed.tmp 1000genomes.line1.chm13v2.0.unmapped
python subscript/mod_label.py 1000genomes.line1.hg38.bed.tmp > 1000genomes.line1.hg38.bed
python subscript/mod_label.py 1000genomes.line1.chm13v2.0.bed.tmp > 1000genomes.line1.chm13v2.0.bed

##########
# gnomAD SV v4.1
[ -f gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz ] || \
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz
[ -f gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz.tbi ] || \
    wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz.tbi

bcftools filter -i 'ALT == "<INS:ME:LINE1>" && INFO/SVLEN >= 5800' gnomad.v4.1.sv.non_neuro_controls.sites.vcf.gz | cut -f 1-8 > gnomad.line1.hg38.vcf
python subscript/proc_gnomad.py gnomad.line1.hg38.vcf > gnomad.line1.hg38.bed

liftOver gnomad.line1.hg38.bed hg38ToHg19.over.chain.gz gnomad.line1.hg19.bed.tmp gnomad.line1.hg19.unmapped
liftOver gnomad.line1.hg38.bed hg38-chm13v2.over.chain.gz gnomad.line1.chm13v2.0.bed.tmp gnomad.line1.chm13v2.0.unmapped
python subscript/mod_label.py gnomad.line1.hg19.bed.tmp > gnomad.line1.hg19.bed
python subscript/mod_label.py gnomad.line1.chm13v2.0.bed.tmp > gnomad.line1.chm13v2.0.bed

##########
# HGSVC3 MEI
python subscript/proc_hgsvc3.py MEI_Callset_GRCh38.ALL.20241211.csv.gz > hgsvc3.line1.hg38.bed

liftOver hgsvc3.line1.hg38.bed hg38ToHg19.over.chain.gz hgsvc3.line1.hg19.bed.tmp hgsvc3.line1.hg19.unmapped
python subscript/mod_label.py hgsvc3.line1.hg19.bed.tmp > hgsvc3.line1.hg19.bed

python subscript/proc_hgsvc3.py MEI_Callset_T2T-CHM13.ALL.20241211.csv.gz > hgsvc3.line1.chm13v2.0.bed

##########
# Merge with priority-based dedup (rmsk > HGSVC3 > 1000genomes > gnomAD)
for ref in hg38 hg19 chm13v2.0; do
    python subscript/merge_bed_with_dedup.py \
        rmsk.line1.${ref}.bed \
        hgsvc3.line1.${ref}.bed \
        1000genomes.line1.${ref}.bed \
        gnomad.line1.${ref}.bed \
        > LINE1.${ref}.bed

    bgzip -c LINE1.${ref}.bed > LINE1.${ref}.bed.gz
    tabix -p bed LINE1.${ref}.bed.gz
done

# mv LINE1.hg38.bed.gz LINE1.hg38.bed.gz.tbi ../nanomonsv/data/
# mv LINE1.hg19.bed.gz LINE1.hg19.bed.gz.tbi ../nanomonsv/data/
# mv LINE1.chm13v2.0.bed.gz LINE1.chm13v2.0.bed.gz.tbi ../nanomonsv/data/

##########
# Cleanup intermediate files
rm -f rmsk.line1.* 1000genomes.line1.* gnomad.line1.* hgsvc3.line1.*
rm -f LINE1.hg38.bed LINE1.hg19.bed LINE1.chm13v2.0.bed
rm -f chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
