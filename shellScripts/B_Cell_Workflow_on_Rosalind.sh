#!/bin/sh
#SBATCH --time=20:00:00
#SBATCH -p shared,brc
#SBATCH --mem 200G

PATH=$R_HOME/bin:$PATH

module load utilities/use.dev
export TMPDIR=/users/k20064105/temp/
module load apps/rstudio/4.1.0-singularity

cd /users/k20064105/flowCytometryDataAnalysis

Rscript R/bCells/01_workflow_preprocessing.R
Rscript R/bCells/02_workflow_convertToDataFrame.R
Rscript R/bCells/03_workflow_multipleRegressionTesting.R
Rscript R/bCells/04_workflow_flowsomClustering.R
Rscript R/bCells/05_workflow_phenographClustering.R
Rscript R/bCells/06_workflow_fastPGClustering.R
#Rscript R/bCells/07_workflow_umapDimReduction.R
#Rscript R/bCells/08_workflow_visuliseUmap.R
#Rscript R/bCells/09_workflow_diffusionMapDimReduction.R
#Rscript R/bCells/10_workflow_visuliseDiffusionMap.R