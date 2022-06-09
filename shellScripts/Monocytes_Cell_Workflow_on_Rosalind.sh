#!/bin/sh
#SBATCH --time=20:00:00
#SBATCH -p shared,brc
#SBATCH --mem 100G

PATH=$R_HOME/bin:$PATH

module load utilities/use.dev
export TMPDIR=/users/k20064105/temp/
module load apps/rstudio/4.1.0-singularity

cd /users/k20064105/flowCytometryDataAnalysis

Rscript R/monocytes/01_workflow_preprocessing.R
Rscript R/monocytes/02_workflow_convertToDataFrame.R
Rscript R/monocytes/04_workflow_flowsomClustering.R