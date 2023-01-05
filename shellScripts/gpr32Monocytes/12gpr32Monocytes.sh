#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 200G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/01_workflow_preprocessing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/02_workflow_convertToDataFrame.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/03_workflow_convertToFcs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/04_workflow_copyToClusteringOutputs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/05_workflow_flowsomClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/06_workflow_flowSomElbowPlot.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/07_workflow_updateFlowSomColumn.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/08_workflow_phenographClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/09_workflow_fastPGClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/10_workflow_identifyCellPopulations.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/11_workflow_umapDimReduction.R
singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/12_workflow_visuliseUmap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/13_workflow_diffusionMapDimReduction.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/14_workflow_visuliseDiffusionMap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/15_calculate_counts.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/16_calculate_unstable_clusters.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/17_flowsom_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/18_reduce_bootstrap_size.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/19_identify_cell_populations_for_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/20_flowsom_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/21_phenograph_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/22_identify_phenograph_cell_populations_for_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/23_phenograph_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/24_differential_abundance_testing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/25_cluster_heatmap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32Monocytes/26_calculate_median_expression.R

sbatch shellScripts/gpr32Monocytes/13gpr32Monocytes.sh