#!/usr/bin/env bash

# SLURM
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=999.logs/003.VCF_QC.log

ml build-env/2020
ml vcftools/0.1.16-foss-2018b-perl-5.28.0

# quality check the VCF file
# based on https://speciationgenomics.github.io/filtering_vcfs/

# DATA #
VCF=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/001.population_structure/001.data/Lactuca.snp.TKI.sub.vcf.gz
OUT=/groups/nordborg/user/pieter.clauw/Documents/Collaborations/salad_marco/001.population_structure/001.data/Lactuca.snp.TKI.sub

#Calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2

Calculate mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT

#Calculate mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT

#Calculate site quality
vcftools --gzvcf $VCF --site-quality --out $OUT

#Calculate proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT

#Calculate proportion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT

#Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $VCF --het --out $OUT



