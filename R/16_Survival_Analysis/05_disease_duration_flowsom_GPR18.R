try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterName <- clusterColumns[3]
markersOrCell <- markersOrCellsClassification[3]
markerName <- "GPR18"

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
    delataAlsfrsRScore - 1"

scaleDatas <- c(TRUE, FALSE)
robusts <- c(TRUE, FALSE)

for (robust in robusts) {
  for (scaleData in scaleDatas) {
    message(paste0("scaleData: ", scaleData))
    message(paste0("robust: ", robust))

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
    message()
  }
}
