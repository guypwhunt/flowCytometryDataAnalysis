#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem=200G

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/01_workflow_preprocessing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/02_workflow_convertToDataFrame.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/03_workflow_convertToFcs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/04_workflow_copyToClusteringOutputs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/05_workflow_flowsomClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/06_workflow_flowSomElbowPlot.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/07_workflow_updateFlowSomColumn.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/08_workflow_phenographClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/09_workflow_fastPGClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/10_workflow_identifyCellPopulations.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/11_workflow_umapDimReduction.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/12_workflow_visuliseUmap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/13_workflow_diffusionMapDimReduction.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/14_workflow_visuliseDiffusionMap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/15_calculate_counts.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/16_calculate_unstable_clusters.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/17_flowsom_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/18_reduce_bootstrap_size.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/19_identify_cell_populations_for_stability_analysis.R
singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/20_flowsom_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/21_phenograph_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/22_identify_phenograph_cell_populations_for_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/23_phenograph_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/24_differential_abundance_testing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/25_cluster_heatmap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr18Senescence/26_calculate_median_expression.R

sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0
sleep 1h
sbatch shellScripts/gpr18Senescence/21gpr18Senescence.sh 0

