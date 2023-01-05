#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH -p singularity_cpu
#SBATCH --mem=1000G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

iteration=${1}
addition=1
limit=11

if [ "$iteration" -lt "$limit" ]; then
  echo "Iteration is less than limit"
  echo $iteration
  singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/senescence/21_phenograph_bootstrapping.R
  sbatch shellScripts/ssenescencePhenographStabilitySingularity.sh $((iteration + addition));
fi