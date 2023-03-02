#!/usr/bin/env bash

# activate conda environment.
conda activate gwas_lactuca

# DATA #
# update to correct working directory
WORKdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/
## general data directory, containing the Lactuca VCF file, phenotypes and ID list
DATAdir=${WORKdir}003.gwas/001.data/
## directory for storing GWAS results
RESULTSdir=${WORKdir}003.gwas/003.results/
## directory for temporary files
TEMP=${WORKdir}003.gwas/temp/

## phenotypes: the first column has the TKI ID for each accession, the second column contains phenotype measurements. Columns are tab separated
PHENO=${DATAdir}leafAnthocyanin.txt
## bed file containing genotype information, created with 001.prepare_bed.sh
BED=${TEMP}Lactuca.snp.TIK.sub

## extract the name of the phenotype (column name)
PHENOname=$(awk '(NR == 1){print $2}' $PHENO)
## temporary file storing the phenotype data
PHENOtemp=${TEMP}${PHENOname}.txt

## temporary file combining genotype and phenotype information. Input for GWAS
GWASinput=${TEMP}$(basename $BED).${PHENOname}
## file storing the kinship matrix, will be saved in GEMMAdir
K=kinship_${PHENOname}
## directory for running GWAS and storing output
GEMMAdir=${RESULTSdir}${PHENOname}_gemma/

# remove header and family column (column of zeroes) to phenoype file
awk '(NR > 1){print 0, $1, $2}' $PHENO > $PHENOtemp

# update .fam file, prepare .bed file for GWAS
plink2 --bfile $BED --pheno $PHENOtemp --make-bed --out $GWASinput

# run PCA that will be used as covariates for population structure correction (in addition to kinship)
## this step takes ±25min.
plink2 --bfile $GWASinput --pca 5 --out $GWASinput
awk '(NR > 1) {print 1, $3, $4, $5, $6, $7}' ${GWASinput}.eigenvec > ${GWASinput}.covariates

# change to GEMMAdir, where output will be written
mkdir -p $GEMMAdir
cd $GEMMAdir

# calculate kinship and save in K
# this can take some time
gemma -bfile $GWASinput -gk 1 -o $K

# run GWAS
gemma -bfile $GWASinput -k output/${K}.cXX.txt -c ${GWASinput}.covariates -lmm 4 -maf 0.05 -o $PHENOname

# cleanup
rm -v ${GWASinput}*
rm -v $PHENOtemp
