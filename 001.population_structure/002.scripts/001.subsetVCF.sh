#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/00x.selectSamplesForManual.log

ml build-env/2020
ml bcftools/1.9-foss-2018b

# subset the SNP matrix for accessions used in our study

# DATA #
VCF=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/Lactuca.snp.vcf.gz
ACCESSIONS=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/Lettuce_accessions_data.txt
IDkey=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/000.generalData/001.data/445LettuceRunningIDs.txt


## replace running IDs by TKI IDs in VCF before subsetting
### create ID map
ID_MAP=/scratch-cbe/users/pieter.clauw/temp/ID_map.txt
VCF_IDs=/scratch-cbe/users/pieter.clauw/temp/VCF_IDs.txt
bcftools query -l $VCF > $VCF_IDs
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' $IDkey $VCF_IDs > $ID_MAP

# get accession IDs for subsetting
IDsub=/scratch-cbe/users/pieter.clauw/temp/ID_sub.txt
awk '(NR > 1) {print $1}' $ACCESSIONS > $IDsub

## update VCF header and subset
bcftools reheader -s $ID_MAP $VCF | bcftools view -Oz -a -S $IDsub -o ${VCF/.vcf.gz/.TKI.sub.vcf.gz}
tabix ${VCF/.vcf.gz/.TKI.sub.vcf.gz}

