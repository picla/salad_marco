#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/002.GWAS_%A_%a.log
#SBATCH --array=3-10

# load required software
ml build-env/2020
ml plink/2.00-alpha1-20190826-avx

# activate GEMMA conda environment
source activate /users/pieter.clauw/.conda/envs/gemma

# DATA #
## this variable will select the phenotype column for GWAS
i=$SLURM_ARRAY_TASK_ID

## general data directory, containing the phenotype data
DATAdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/
## directory for storing GWAS results
RESULTSdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/002.GWAS_test/003.results/
## working directory
SCRATCHdir=/scratch-cbe/users/pieter.clauw/salad_marco/002.GWAS_test/

## phenotypes with first column containing the TKI ID for each accession, subsequent column(s) contain phenotype measurements
PHENO=${DATAdir}tableS10_agronomicTraits.txt
## bed file containing genotype information, created with 001.prepare_bed.sh
BED=${SCRATCHdir}001.data/Lactuca.snp.TIK.sub
## directory for temporary files
TEMP=/scratch-cbe/users/pieter.clauw/temp/


## extract the name of the phenotype (column name)
PHENOname=$(awk -v i=$i 'BEGIN {FS = "\t"}(NR == 1){print $i}' $PHENO | sed "s/ /_/g")
## temporary file storing the phenptype data
PHENOtemp=${TEMP}${PHENOname}.txt

## temporary file combining genotype and phenotype information. Input for GWAS
GWASinput=${TEMP}$(basename $BED).${PHENOname}
## file storing the kinship matrix, will be saved in GEMMAdir
K=kinship_${PHENOname}
## directory for running GWAS and storing output
GEMMAdir=${SCRATCHdir}003.results/GEMMA/

## output whihc phenotype is being used
echo job $i corresponds to phenotype $PHENOname

# select phenotype and store in temporary file
awk -v i=$i '(NR > 1){print 0, $1, $i}' $PHENO > $PHENOtemp
## remove NA values
sed -i 's/ -/ NA/g' $PHENOtemp

# update .fam file, prepare .bed file for GWAS
plink2 --bfile $BED --pheno $PHENOtemp --make-bed --out $GWASinput

# run PCA that will be used as covariates for population structure correction (in addition to kinship)
plink2 --bfile $GWASinput --pca 5 --out $GWASinput
awk '(NR > 1) {print 1, $3, $4, $5, $6, $7}' ${GWASinput}.eigenvec > ${GWASinput}.covariates

# change to GEMMAdir, where output will be written
mkdir -p $GEMMAdir
cd $GEMMAdir

# calculate kinship and save in K
gemma -bfile $GWASinput -gk 1 -o $K

# run GWAS
gemma -bfile $GWASinput -k output/${K}.cXX.txt -c ${GWASinput}.covariates -lmm 4 -maf 0.05 -o $PHENOname

# cleanup
rm -v ${GWASinput}*
rm -v $PHENOtemp

# copy to results to the RESULTSdir
cp ${GEMMAdir}output/${PHENOname}* $RESULTSdir