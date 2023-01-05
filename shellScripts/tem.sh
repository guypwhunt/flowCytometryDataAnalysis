#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH -p cpu
#SBATCH --mem 20G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0
sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0
sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0
sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0
sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0
sleep 1h
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Monocytes/21gpr32Monocytes.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32Senescence/21gpr32Senescence.sh 0
sbatch /scratch/users/k20064105/flowCytometryDataAnalysis/shellScripts/gpr32TCells/21gpr32TCells.sh 0