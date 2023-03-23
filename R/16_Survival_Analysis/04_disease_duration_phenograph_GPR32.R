try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterName <- clusterColumns[4]
markersOrCell <- markersOrCellsClassification[3]
markerName <- "GPR32"

fileName <- paste0(markerName,
                   clusterName,
                   "DiseaseDurationBiologicalModel")

clinicalCovariates <- c("ethnicityID",
                        "alsfrsR",
                        "ageAtOnset",
                        "riluzole",
                        "delataAlsfrsRScore")

clinicalColumnsToScale <- c("alsfrsR",
                            "ageAtOnset",
                            "delataAlsfrsRScore")

coxModelForumla <-
  "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
    riluzole + delataAlsfrsRScore + GPR32_median_Positive_Follicular_B_Cells +
    GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_. +
    GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytes +
    GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low. +
    GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytes +
    GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low. -
    1"

scaleDatas <- c(TRUE, FALSE)
robusts <- c(TRUE, FALSE)

for (robust in robusts) {
  for (scaleData in scaleDatas) {
    minDF <-
      survivalModelAnalysis(
        fileName,
        clinicalCovariates,
        clinicalColumnsToScale,
        scaleData = scaleData,
        robust = robust,
        clusterName = clusterName,
        markersOrCell = markersOrCell,
        markerName = markerName
      )

    saveSurvivalModel(minDF,
                      coxModelForumla,
                      fileName,
                      scaleData = scaleData,
                      robust = robust)
  }
}
