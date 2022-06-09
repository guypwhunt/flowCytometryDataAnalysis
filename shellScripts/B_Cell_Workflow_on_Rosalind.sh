#!/bin/sh
#SBATCH --time=02:00:00
#SBATCH -p shared,brc
#SBATCH --mem 100G

PATH=$R_HOME/bin:$PATH

module load utilities/use.dev
export TMPDIR=/users/k20064105/temp/
module load apps/rstudio/4.1.0-singularity

cd /users/k20064105/flowCytometryDataAnalysis

Rscript R/bCells/04_workflow_flowsomClustering.R
Rscript R/bCells/11_differential_abundance_testing_test.R