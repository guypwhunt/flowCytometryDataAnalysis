try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterNames <- clusterColumns[3:4]

markersOrCell <- markersOrCellsClassification[3]
markerNames <- c("GPR18", "GPR32")

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo <- experimentInfo[experimentInfo$visit == 1,]
experimentInfo <-
  experimentInfo[experimentInfo$experiment == "flowCytometry",]
experimentInfo <-
  experimentInfo[experimentInfo$caseControl == "Case",]

filePath <- "data/medianValues/"

for(clusterName in clusterNames) {
  fileNames <- list.files(filePath)
  fileNames <-
    fileNames[grep(markersOrCell, fileNames, fixed = FALSE)]
  fileNames <-
    fileNames[grep(clusterName, fileNames, fixed = FALSE)]

  fileNamesPath <- paste0(filePath, fileNames)

  dfs <- lapply(fileNamesPath, read.csv)

  for (df in dfs) {
    colnames(df)[2:ncol(df)] <- paste0(colnames(df)[2:ncol(df)], clusterName)
    if (exists("mergedDF")) {
      mergedDF <- merge(mergedDF, df, by = "fileName", all = TRUE)
    } else {
      mergedDF <- df
    }
  }
}

mergedDF <- mergedDF[mergedDF$fileName %in% experimentInfo$sample_id,]

message(nrow(mergedDF))
mergedDF  <- mergedDF[, colSums(is.na(mergedDF))<ceiling(nrow(mergedDF)/10)]
message(nrow(mergedDF))

combinedDf <- merge(
  mergedDF,
  experimentInfo,
  by.x = "fileName",
  by.y = "sample_id",
  all.x = TRUE
)

dateOfLateOutCome <-
  combinedDf[is.na(combinedDf$outcomeDeathDate), c("outcomeLastVisitDate", "visitDate")]

dateOfLateOutCome <-
  (dateOfLateOutCome$outcomeLastVisitDate - dateOfLateOutCome$visitDate) /
  (24 * 60 * 60 * 365.25)

combinedDf[is.na(combinedDf$outcomeDeathDate), "diseaseDurationInYears"] <-
  dateOfLateOutCome

combinedDf$fileName <- factor(combinedDf$fileName)

combinedDf$diseaseDurationInYears <-
  combinedDf$diseaseDurationInYears * 365.25

combinedDf$outcome <- factor(combinedDf$outcome)
combinedDf$status <- as.numeric(combinedDf$outcome)
combinedDf$status <- ifelse(combinedDf$status == 2, 1, 0)

combinedDf$gender <- factor(combinedDf$gender)
combinedDf$sex <- as.numeric(combinedDf$gender)

combinedDf$ethnicity <- factor(combinedDf$ethnicity)
combinedDf$ethnicityID <- as.numeric(combinedDf$ethnicity)

combinedDf$BulbarLimb <- factor(combinedDf$BulbarLimb)
combinedDf$onset <- as.numeric(combinedDf$BulbarLimb)

clinicalCovariates <- c("ethnicityID",
                        "alsfrsR",
                        "ageAtOnset",
                        "riluzole",
                        "delataAlsfrsRScore")

biologicalCovariates <- c(colnames(mergedDF)[colnames(mergedDF) != "fileName"])

clinicalAndBiologicalCovariates <-
  append(clinicalCovariates, biologicalCovariates)

combinedDf$sex <- as.factor(combinedDf$sex)
combinedDf$ethnicityID <- as.factor(combinedDf$ethnicityID)
combinedDf$onset <- as.factor(combinedDf$onset)

# Censored Analysis
minDF <-
  combinedDf[, c("diseaseDurationInYears",
                 "status",
                 clinicalAndBiologicalCovariates)]
minDF <- na.omit(minDF)


clinical <- readRDS("data/survivalAnalysis/DiseaseDurationClinicalModel.rds")
clinical <-
  coxph(as.formula(clinical$formula)
  , data = minDF)


flowsomGpr32 <- readRDS("data/survivalAnalysis/GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.rds")
print(flowsomGpr32$formula)
flowsomGpr32 <-
  coxph(as.formula(Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     ageAtOnset + delataAlsfrsRScore + GPR32_median_Positive_Follicular_B_Cellsmeta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.meta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.meta_clusters_flowsom -
                     1)
        , data = minDF)

phenographGpr32 <- readRDS("data/survivalAnalysis/GPR32clusters_phenographDiseaseDurationBiologicalModel.rds")
print(phenographGpr32$formula)

phenographGpr32 <-
  coxph(as.formula(Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     riluzole + delataAlsfrsRScore + GPR32_median_Positive_Follicular_B_Cellsclusters_phenograph +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.clusters_phenograph -
                     1)
        , data = minDF)

flowsomGpr18 <- readRDS("data/survivalAnalysis/GPR18meta_clusters_flowsomDiseaseDurationBiologicalModel.rds")
flowsomGpr18 <-
  coxph(as.formula(flowsomGpr18$formula)
        , data = minDF)

phenographGpr18 <- readRDS("data/survivalAnalysis/GPR18clusters_phenographDiseaseDurationBiologicalModel.rds")
phenographGpr18 <-
  coxph(as.formula(phenographGpr18$formula)
        , data = minDF)

flowsomGpr18Gpr32 <- readRDS("data/survivalAnalysis/GPR18GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.rds")
print(flowsomGpr18Gpr32$formula)
flowsomGpr18Gpr32 <-
  coxph(as.formula(Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     ageAtOnset + delataAlsfrsRScore + GPR32_median_Positive_Naïve_B_Cellsmeta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.meta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.meta_clusters_flowsom -
                     1)
        , data = minDF)

phenographGpr18Gpr32 <- readRDS("data/survivalAnalysis/GPR18GPR32clusters_phenographDiseaseDurationBiologicalModel.rds")
print(phenographGpr18Gpr32$formula)
phenographGpr18Gpr32 <-
  coxph(as.formula(Surv(diseaseDurationInYears, status) ~ ethnicityID + ageAtOnset +
                     riluzole + delataAlsfrsRScore + GPR18_median_Positive_Follicular_B_Cellsclusters_phenograph +
                     GPR18_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.clusters_phenograph +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.clusters_phenograph -
                     1)
        , data = minDF)

for (model in list(flowsomGpr18, phenographGpr18, flowsomGpr32, phenographGpr32, flowsomGpr18Gpr32, phenographGpr18Gpr32))  {
  if (exists("ddComparison")) {
    ddComparison <- rbind(ddComparison, as.data.frame(anova(clinical, model))[2,])
  } else {
    ddComparison <- as.data.frame(anova(clinical, model))
  }
}

# ddComparison <-
#   anova(clinical, flowsomGpr18, phenographGpr18, flowsomGpr32, phenographGpr32, flowsomGpr18Gpr32, phenographGpr18Gpr32)

dir.create("data/survivalAnalysis/compareAllModels/")

rownames(ddComparison) <- seq(nrow(ddComparison))

write.csv(
  as.data.frame(ddComparison),
  paste0(
    "data/survivalAnalysis/compareAllModels/",
    "diseaseDuration.csv"
  ),
  row.names = TRUE
)
