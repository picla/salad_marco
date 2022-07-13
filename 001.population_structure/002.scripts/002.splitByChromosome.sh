#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/002.selectSamplesForManual.log
#SBATCH --array=1-9

ml build-env/2020
ml bcftools/1.9-foss-2018b

VCF=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/001.population_structure/001.data/Lactuca.snp.TKI.sub.vcf.gz
chr=$SLURM_ARRAY_TASK_ID
VCFchr=${VCF/.vcf.gz/.chr${chr}.vcf.gz}

bcftools filter -r chr${chr} -Oz -o $VCFchr $VCF
tabix $VCFchr