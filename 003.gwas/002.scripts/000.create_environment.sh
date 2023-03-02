# create a conda environment
# this environment contains GEMMA and PLINK2
# these two software packages are needed to run GWAS
# in case conda is not yet available, check link below for installation guides
# https://conda.io/projects/conda/en/latest/user-guide/install/index.html

git clone https://github.com/picla/salad_marco/
cd salad_marco
conda env create -f environment.yml

