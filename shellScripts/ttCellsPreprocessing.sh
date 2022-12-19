#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 200G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/tCells/24_differential_abundance_testing.R