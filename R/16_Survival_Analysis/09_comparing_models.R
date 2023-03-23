try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterNames <- clusterColumns[3:4]

markersOrCell <- markersOrCellsClassification[3]
markerNames <- c("GPR18", "GPR32")

folderNames <- c("SurvivalAnalysis", "robustSurvivalAnalysis",
               "ScaledSurvivalAnalysis", "robustScaledSurvivalAnalysis")

clinicalColumnsToScale <- c("alsfrsR",
                            "ageAtOnset",
                            "delataAlsfrsRScore")

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx") %>%
  as.data.frame() %>%
  filter(visit == 1) %>%
  filter(experiment == "flowCytometry") %>%
  filter(caseControl == "Case")

filePath <- "data/medianValues/"

clinicalCovariates <- c("ethnicityID",
                        "alsfrsR",
                        "ageAtOnset",
                        "riluzole",
                        "delataAlsfrsRScore")

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

mergedDF  <- mergedDF[, colSums(is.na(mergedDF))<ceiling(nrow(mergedDF)/10)]

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

goldenSource <- minDF

for (folderName in folderNames) {
  try({
    minDF <- goldenSource

    folderPath <- paste0("data/", folderName)

    robust <- grepl("robust", folderName)

    scaleData <- grepl("Scaled", folderName)

    if (scaleData) {
      minDF[,c(clinicalColumnsToScale, biologicalCovariates)] <-
        scale(minDF[,c(clinicalColumnsToScale, biologicalCovariates)])
    }


    clinical <-
      readRDS(paste0("data/", folderName, "/DiseaseDurationClinicalModel.rds"))
    clinical <-
      coxph(as.formula(clinical$formula)
            ,
            data = minDF,
            robust = robust)


    flowsomGpr32 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.rds"
        )
      )
    flowsomGpr32 <-
      coxph(
        as.formula(
          "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     ageAtOnset + delataAlsfrsRScore + GPR32_median_Positive_Follicular_B_Cellsmeta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.meta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.meta_clusters_flowsom -
                     1"
        )
        ,
        data = minDF,
        robust = robust
      )

    phenographGpr32 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR32clusters_phenographDiseaseDurationBiologicalModel.rds"
        )
      )
    phenographGpr32 <-
      coxph(
        as.formula(
          "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     riluzole + delataAlsfrsRScore + GPR32_median_Positive_Follicular_B_Cellsclusters_phenograph +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.clusters_phenograph -
                     1"
        )
        ,
        data = minDF,
        robust = robust
      )

    flowsomGpr18 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR18meta_clusters_flowsomDiseaseDurationBiologicalModel.rds"
        )
      )
    flowsomGpr18 <-
      coxph(as.formula(flowsomGpr18$formula)
            ,
            data = minDF,
            robust = robust)

    phenographGpr18 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR18clusters_phenographDiseaseDurationBiologicalModel.rds"
        )
      )
    phenographGpr18 <-
      coxph(as.formula(phenographGpr18$formula)
            ,
            data = minDF,
            robust = robust)

    flowsomGpr18Gpr32 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR18GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.rds"
        )
      )
    flowsomGpr18Gpr32 <-
      coxph(
        as.formula(
          "Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
                     ageAtOnset + delataAlsfrsRScore + GPR32_median_Positive_Naïve_B_Cellsmeta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.meta_clusters_flowsom +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.meta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesmeta_clusters_flowsom +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.meta_clusters_flowsom -
                     1"
        )
        ,
        data = minDF,
        robust = robust
      )

    phenographGpr18Gpr32 <-
      readRDS(
        paste0(
          "data/",
          folderName,
          "/GPR18GPR32clusters_phenographDiseaseDurationBiologicalModel.rds"
        )
      )
    phenographGpr18Gpr32 <-
      coxph(
        as.formula(
          "Surv(diseaseDurationInYears, status) ~ ethnicityID + ageAtOnset +
                     riluzole + delataAlsfrsRScore + GPR18_median_Positive_Follicular_B_Cellsclusters_phenograph +
                     GPR18_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.clusters_phenograph +
                     GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.clusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytesclusters_phenograph +
                     GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.clusters_phenograph -
                     1"
        )
        ,
        data = minDF,
        robust = robust
      )

    for (model in list(
      flowsomGpr18,
      phenographGpr18,
      flowsomGpr32,
      phenographGpr32,
      flowsomGpr18Gpr32,
      phenographGpr18Gpr32
    ))  {
      if (exists("ddComparison")) {
        ddComparison <-
          rbind(ddComparison, as.data.frame(anova(clinical, model))[2, ])
      } else {
        ddComparison <- as.data.frame(anova(clinical, model))
      }
    }

    dir.create(paste0("data/", folderName, "/compareAllModels/"),
               showWarnings = FALSE)

    rownames(ddComparison) <- seq(nrow(ddComparison))

    write.csv(
      as.data.frame(ddComparison),
      paste0(
        "data/",
        folderName,
        "/compareAllModels/",
        "diseaseDuration.csv"
      ),
      row.names = TRUE
    )

    rm(ddComparison)
  })
}
