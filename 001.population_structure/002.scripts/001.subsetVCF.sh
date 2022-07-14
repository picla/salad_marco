#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/001.subsetVCF.log

ml build-env/2020
ml bcftools/1.9-foss-2018b

# subset the SNP matrix for accessions used in our study

# DATA #
DIRdata=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/
DIRout=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/001.population_structure/001.data/
TEMP=/scratch-cbe/users/pieter.clauw/temp/
VCF=Lactuca.snp.vcf.gz
VCFout=${VCF/.vcf.gz/.TKI.sub.vcf.gz}
ACCESSIONS=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/Lettuce_accessions_data.txt
IDkey=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/445LettuceRunningIDs.txt


## replace running IDs by TKI IDs in VCF before subsetting
### create ID map
ID_MAP=${TEMP}ID_map.txt
VCF_IDs=${TEMP}VCF_IDs.txt
bcftools query -l ${DIRdata}${VCF} > $VCF_IDs
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' $IDkey $VCF_IDs > $ID_MAP

# get accession IDs for subsetting
IDsub=${TEMP}ID_sub.txt
awk '(NR > 1) {print $1}' $ACCESSIONS > $IDsub

## update VCF header and subset
bcftools reheader -s $ID_MAP ${DIRdata}${VCF} | bcftools view -Ou -a -S $IDsub | bcftools view -Oz -m 2 -o ${DIRout}${VCFout}
tabix ${DIRout}${VCFout}

