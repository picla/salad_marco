#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/001.prepare_bed.log

ml build-env/2020
ml plink/2.00-alpha1-20190826-avx
ml bcftools/1.9-foss-2018b

# DATA #
DATAdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/
SCRATCHdir=/scratch-cbe/users/pieter.clauw/salad_marco/002.GWAS_test/
TEMP=/scratch-cbe/users/pieter.clauw/temp/

VCF=${DATAdir}Lactuca.snp.vcf.gz
PHENO=${DATAdir}tableS10_agronomicTraits.txt
IDkey=${DATAdir}445LettuceRunningIDs.txt

ID_MAP=${TEMP}ID_map.txt
VCF_IDs=${TEMP}VCF_IDs.txt
IDsub=${TEMP}IDsub.txt

VCFsub=${SCRATCHdir}001.data/$(basename -s .vcf.gz $VCF).TIK.sub.vcf.gz

# filter VCF for phenotyped accessions
## get VCF sample names
bcftools query -l $VCF > $VCF_IDs
## link VCF sample names to TIK IDs
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' $IDkey $VCF_IDs > $ID_MAP
# get phenotyped accessions
awk '(NR > 1) {print $1}' $PHENO > $IDsub
## rename VCF sample names according ID_MAP and subset for phenotyped accessions
bcftools reheader -s $ID_MAP $VCF | bcftools view -a -S $IDsub -Oz -o $VCFsub
tabix $VCFsub

# create .bed .fam and .bim files
plink2 --vcf $VCFsub --make-bed --out ${VCFsub/.vcf.gz/}

# cleanup
rm -v $ID_MAP
rm -v $VCF_IDs
rm -v $IDsub
