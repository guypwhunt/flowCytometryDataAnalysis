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

coxModelForumla <-
  "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
        ageAtOnset + riluzole + delataAlsfrsRScore - 1"

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
        robust = robust
      )

    saveSurvivalModel(minDF,
                      coxModelForumla,
                      fileName,
                      scaleData = scaleData,
                      robust = robust)
  }
}
