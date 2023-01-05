#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH -p cpu
#SBATCH --mem=200G

PATH=$R_HOME/bin:$PATH

cd /scratch/users/k20064105/flowCytometryDataAnalysis

##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/01_workflow_preprocessing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/02_workflow_convertToDataFrame.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/03_workflow_convertToFcs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/04_workflow_copyToClusteringOutputs.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/05_workflow_flowsomClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/06_workflow_flowSomElbowPlot.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/07_workflow_updateFlowSomColumn.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/08_workflow_phenographClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/09_workflow_fastPGClustering.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/10_workflow_identifyCellPopulations.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/11_workflow_umapDimReduction.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/12_workflow_visuliseUmap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/13_workflow_diffusionMapDimReduction.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/14_workflow_visuliseDiffusionMap.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/15_calculate_counts.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/16_calculate_unstable_clusters.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/17_flowsom_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/18_reduce_bootstrap_size.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/19_identify_cell_populations_for_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/20_flowsom_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/21_phenograph_bootstrapping.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/22_identify_phenograph_cell_populations_for_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/23_phenograph_cluster_stability_analysis.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/24_differential_abundance_testing.R
##singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/25_cluster_heatmap.R
singularity exec --home /scratch/users/k20064105/flowCytometryDataAnalysis rstudio_4-2-0_latest.sif Rscript R/gpr32TCells/26_calculate_median_expression.R