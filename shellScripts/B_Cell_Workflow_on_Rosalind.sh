#!/bin/sh
#SBATCH --time=20:00:00
#SBATCH -p shared,brc
#SBATCH --mem 300G

PATH=$R_HOME/bin:$PATH

module load utilities/use.dev
export TMPDIR=/users/k20064105/temp/
module load apps/rstudio/4.1.0-singularity

cd /users/k20064105/flowCytometryDataAnalysis

Rscript R/bcells/01_workflow_preprocessing.R
Rscript R/bcells/02_workflow_convertToDataFrame.R
Rscript R/bcells/03_workflow_multipleRegressionTesting.R
Rscript R/bcells/04_workflow_flowsomClustering.R
Rscript R/bcells/05_workflow_phenographClustering.R
Rscript R/bcells/06_workflow_fastPGClustering.R
Rscript R/bcells/07_workflow_umapDimReduction.R
Rscript R/bcells/08_workflow_visuliseUmap.R
Rscript R/bcells/09_workflow_diffusionMapDimReduction.R
Rscript R/bcells/10_workflow_visuliseDiffusionMap.R