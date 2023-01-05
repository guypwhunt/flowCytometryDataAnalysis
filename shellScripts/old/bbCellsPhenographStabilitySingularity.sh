#!/bin/sh
#SBATCH --time=36:00:00
#SBATCH -p singularity_cpu
#SBATCH --mem 300G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

iteration=${1}
addition=1
limit=3

if [ "$iteration" -lt "$limit" ]; then
  echo "Iteration is less than limit"
  echo $iteration
  singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/bCells/21_phenograph_bootstrapping.R
  sbatch shellScripts/bbCellsPhenographStabilitySingularity.sh $((iteration + addition));
fi