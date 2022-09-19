#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/002.GWAS_%A_%a.log
#SBATCH --array=3-10

ml build-env/2020
ml plink/2.00-alpha1-20190826-avx

# activate GEMMA conda environment
source activate /users/pieter.clauw/.conda/envs/gemma

# DATA #
i=$SLURM_ARRAY_TASK_ID

DATAdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/
RESULTSdir=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/002.GWAS_test/003.results/
SCRATCHdir=/scratch-cbe/users/pieter.clauw/salad_marco/002.GWAS_test/

PHENO=${DATAdir}tableS10_agronomicTraits.txt
BED=${SCRATCHdir}001.data/Lactuca.snp.TIK.sub
TEMP=/scratch-cbe/users/pieter.clauw/temp/

PHENOname=$(awk -v i=$i 'BEGIN {FS = "\t"}(NR == 1){print $i}' $PHENO | sed "s/ /_/g")
PHENOtemp=${TEMP}${PHENOname}.txt

GWASinput=${TEMP}$(basename $BED).${PHENOname}
K=kinship_${PHENOname}
GEMMAdir=${SCRATCHdir}003.results/GEMMA/

echo job $i corresponds to phenotype $PHENOname

# select phenotype
awk -v i=$i '(NR > 1){print 0, $1, $i}' $PHENO > $PHENOtemp
sed -i 's/ -/ NA/g' $PHENOtemp

# update .fam file, prepare .bed file for GWAS
plink2 --bfile $BED --pheno $PHENOtemp --make-bed --out $GWASinput

# PCA for covariates
plink2 --bfile $GWASinput --pca 5 --out $GWASinput

awk '(NR > 1) {print 1, $3, $4, $5, $6, $7}' ${GWASinput}.eigenvec > ${GWASinput}.covariates

# change to GEMMAdir, where output will be written
mkdir -p $GEMMAdir
cd $GEMMAdir

# kinship
gemma -bfile $GWASinput -gk 1 -o $K

# GWAS
gemma -bfile $GWASinput -k output/${K}.cXX.txt -c ${GWASinput}.covariates -lmm 4 -maf 0.05 -o $PHENOname

# cleanup
rm -v ${GWASinput}*
rm -v $PHENOtemp

# copy to projects
cp ${GEMMAdir}output/${PHENOname}* $RESULTSdir