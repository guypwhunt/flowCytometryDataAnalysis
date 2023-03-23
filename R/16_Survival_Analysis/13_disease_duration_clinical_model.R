try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

fileName <- "DiseaseDurationClinicalModel"

clinicalCovariates <- c(
  "sex",
  "ethnicityID",
  "onset",
  "alsfrsR",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "riluzole",
  "delataAlsfrsRScore"
)

clinicalColumnsToScale <- c("alsfrsR",
                    "diagnosticDelayInYears",
                    "ageAtOnset",
                    "delataAlsfrsRScore")

coxModelForumla <- "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
        ageAtOnset + riluzole + delataAlsfrsRScore - 1"

minDF <- survivalModelAnalysis(fileName, clinicalCovariates, clinicalColumnsToScale, scaleData = FALSE)

saveSurvivalModel(minDF, coxModelForumla, fileName, scaleData = FALSE)


##############
clusterName <- clusterColumns[3]
markersOrCell <- markersOrCellsClassification[3]
markerName <- "GPR32"

minDF <- survivalModelAnalysis(fileName, clinicalCovariates, clinicalColumnsToScale, scaleData = TRUE, clusterName, markersOrCell, markerName)

saveSurvivalModel(minDF, coxModelForumla, fileName, scaleData = TRUE)
