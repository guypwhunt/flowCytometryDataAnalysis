#!/bin/sh
#SBATCH --time=03:00:00
#SBATCH -p shared,brc
#SBATCH --mem 100G

PATH=$R_HOME/bin:$PATH

module load utilities/use.dev
export TMPDIR=/users/k20064105/temp/
module load apps/rstudio/4.1.0-singularity

cd /users/k20064105/flowCytometryDataAnalysis

Rscript R/senescence/04_workflow_flowsomClustering.R