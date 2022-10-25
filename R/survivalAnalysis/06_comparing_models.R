try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

cpddbm <-
  readRDS("data/survivalAnalysis/clusters_phenographDiseaseDurationBiologicalModel.rds")
cpddcm <-
  readRDS("data/survivalAnalysis/clusters_phenographDiseaseDurationClinicalModel.rds")
cpsfvbm <-
  readRDS(
    "data/survivalAnalysis/clusters_phenographSurvivalFromVisitBiologicalModel.rds"
  )
cpsfvcm <-
  readRDS("data/survivalAnalysis/clusters_phenographSurvivalFromVisitClinicalModel.rds")
mcfddbm <-
  readRDS(
    "data/survivalAnalysis/meta_clusters_flowsomDiseaseDurationBiologicalModel.rds"
  )
mcfddcm <-
  readRDS("data/survivalAnalysis/meta_clusters_flowsomDiseaseDurationClinicalModel.rds")
mcfsfvbm <-
  readRDS(
    "data/survivalAnalysis/meta_clusters_flowsomSurvivalFromVisitBiologicalModel.rds"
  )
mcfsfvcm <-
  readRDS(
    "data/survivalAnalysis/meta_clusters_flowsomSurvivalFromVisitClinicalModel.rds"
  )

ddComparison <-
  anova(cpddcm, mcfddcm, cpddbm, mcfddbm
  )

dir.create("data/survivalAnalysis/compareAllModels/")

write.csv(
  as.data.frame(ddComparison),
  paste0(
    "data/survivalAnalysis/compareAllModels/",
    "diseaseDuration.csv"
  ),
  row.names = TRUE
)

sfvComparison <-
  anova(cpsfvcm, mcfsfvcm, cpsfvbm, mcfsfvbm
  )

dir.create("data/survivalAnalysis/compareAllModels/")

write.csv(
  as.data.frame(sfvComparison),
  paste0(
    "data/survivalAnalysis/compareAllModels/",
    "survivalFromVisit.csv"
  ),
  row.names = TRUE
)
