#!/usr/bin/env bash

# activate conda environment.
conda activate gwas_lactuca

# DATA #
## general data directory, containing the Lactuca VCF file, phenotypes and ID list
DATAdir=003.gwas/001.data/
## directory for temporary files
TEMP=003.gwas/temp/

## phenotypes: the first column has the TKI ID for each accession, the second column contains phenotype measurements. Columns are tab separated
PHENO=${DATAdir}leafAnthocyanin.txt
## The general Lactuca VCF file (https://ftp.cngb.org/pub/CNSA/data2/CNP0000335/Other/variation/Lactuca.snp.vcf.gz)
VCF=${DATAdir}Lactuca.snp.vcf.gz
## annotation file for running IDs and TKI IDs
IDkey=${DATAdir}445LettuceRunningIDs.txt

## temporary file storing VCF sample names
VCF_IDs=${TEMP}VCF_IDs.txt
## temporary file to link VCF sample names to TKI IDs
ID_MAP=${TEMP}ID_map.txt
## temporary file with the list of phenotyped accessions
IDsub=${TEMP}IDsub.txt
## subset of VCF file containing only phenotyped accessions
VCFsub=${TEMP}$(basename -s .vcf.gz $VCF).TIK.sub.vcf.gz

# filter VCF for phenotyped accessions
## get VCF sample names
bcftools query -l $VCF > $VCF_IDs
## link VCF sample names to TKI IDs
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' $IDkey $VCF_IDs > $ID_MAP
# get phenotyped accessions
awk '(NR > 1) {print $1}' $PHENO > $IDsub
## rename VCF sample names according ID_MAP and subset for phenotyped accessions
## this takes time (±1 hour)
bcftools reheader -s $ID_MAP $VCF | bcftools view -a -S $IDsub -Oz -o $VCFsub
tabix $VCFsub

# create .bed .fam and .bim files, required for GWAS.
plink2 --vcf $VCFsub --make-bed --out ${VCFsub/.vcf.gz/}

# cleanup
rm -v $ID_MAP
rm -v $VCF_IDs
rm -v $IDsub
