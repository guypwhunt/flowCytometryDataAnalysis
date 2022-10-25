try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification

clusterName <- clusterNames[1]
markersOrCell <- markersOrCells[1]

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo <- experimentInfo[experimentInfo$visit == 1,]
experimentInfo <-
  experimentInfo[experimentInfo$experiment == "flowCytometry",]
experimentInfo <-
  experimentInfo[experimentInfo$caseControl == "Case",]

filePath <- "data/medianValues/"

fileNames <- list.files(filePath)
fileNames <-
  fileNames[grep(clusterName, fileNames, fixed = FALSE)]
fileNames <-
  fileNames[grep(markersOrCell, fileNames, fixed = FALSE)]

fileNamesPath <- paste0(filePath, fileNames)

dfs <- lapply(fileNamesPath, read.csv)

for (df in dfs) {
  if (exists("mergedDF")) {
    mergedDF <- merge(mergedDF, df, by = "fileName", all = TRUE)
  } else {
    mergedDF <- df
  }
}

#mergedDF <- mergedDF[mergedDF$fileName %in% experimentInfo$sample_id,]
#experimentInfo <- experimentInfo[experimentInfo$sample_id %in% mergedDF$fileName,]

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

clinicalCovariates <- c(
  "ageAtVisit",
  "sex",
  "ethnicityID",
  "onset",
  "alsfrsR",
  "timeFromOnsetToVisitInYears",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "riluzole",
  "delataAlsfrsRScore"
)

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

# Univariate Analysis
univ_formulas <- sapply(clinicalAndBiologicalCovariates,
                        function(x)
                          as.formula(paste(
                            'Surv(diseaseDurationInYears, status)~', x
                          )))

univ_models <-
  lapply(univ_formulas, function(x) {
    coxph(x, data = minDF)
  })
# Extract data
univ_results <- lapply(univ_models,
                       function(x) {
                         x <- summary(x)
                         p.value <-
                           signif(x$wald["pvalue"], digits = 2)
                         wald.test <-
                           signif(x$wald["test"], digits = 2)
                         beta <-
                           signif(x$coef[1], digits = 2)
                         #coeficient beta
                         HR <-
                           signif(x$coef[2], digits = 2)
                         #exp(beta)
                         HR.confint.lower <-
                           signif(x$conf.int[, "lower .95"], 2)
                         HR.confint.upper <-
                           signif(x$conf.int[, "upper .95"], 2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res <- c(beta, HR, wald.test, p.value)
                         names(res) <-
                           c("beta", "HR (95% CI for HR)", "wald.test",
                             "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

univ_results$ethnicityID <- univ_results$ethnicityID[c(1, 2, 4, 5)]
names(univ_results$ethnicityID) <- names(univ_results$sex)
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

res.cox <-
  coxph(as.formula(
    paste(
      "Surv(diseaseDurationInYears, status) ~ 0",
      sapply(list(clinicalCovariates), paste, collapse =
               " + "),

      sep = " + "
    )
  )
  , data = minDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

censored.clinical.res.cox <-
  coxph(
    Surv(diseaseDurationInYears, status) ~ alsfrsR + timeFromOnsetToVisitInYears -
      1,
    data = minDF
  )
censored.clinical.res.cox.summary <-
  summary(censored.clinical.res.cox)
saveRDS(censored.clinical.res.cox, file = paste0(
  "data/survivalAnalysis/",
  clusterName,
  "DiseaseDurationClinicalModel.rds"
))
write.csv(
  as.data.frame(censored.clinical.res.cox.summary$coefficients),
  paste0(
    "data/survivalAnalysis/",
    clusterName,
    "DiseaseDurationClinicalModel.csv"
  ),
  row.names = TRUE
)

ggsurvplot(
  survfit(censored.clinical.res.cox, data = minDF),
  palette = "#2E9FDF",
  ggtheme = theme_minimal()
)

res.cox <-
  coxph(as.formula(
    paste(
      "Surv(diseaseDurationInYears, status) ~ 0",
      sapply(list(clinicalAndBiologicalCovariates), paste, collapse =
               " + "),

      sep = " + "
    )
  )
  , data = minDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

censored.clinical.biological.res.cox <-
  coxph(
    Surv(diseaseDurationInYears, status) ~ sex + onset + alsfrsR +
      timeFromOnsetToVisitInYears + diagnosticDelayInYears + ageAtOnset +
      riluzole + delataAlsfrsRScore + GPR32_median_Positive_Naïve_B_Cells +
      GPR32_median_Positive_HLA_Negative_DR_Negative_Or_Low_Activated_CD11b_Positive__Classical_Monocytes +
      GPR32_median_Positive_HLA_Negative_DR_Negative_Or_Low_Intermediate_Monocytes -
      1,
    data = minDF
  )
saveRDS(censored.clinical.biological.res.cox, file = paste0(
  "data/survivalAnalysis/",
  clusterName,
  "DiseaseDurationBiologicalModel.rds"
))
censored.clinical.biological.res.cox.summary <-
  summary(censored.clinical.biological.res.cox)
write.csv(
  as.data.frame(censored.clinical.biological.res.cox.summary$coefficients),
  paste0(
    "data/survivalAnalysis/",
    clusterName,
    "DiseaseDurationBiologicalModel.csv"
  ),
  row.names = TRUE
)

ggsurvplot(
  survfit(censored.clinical.biological.res.cox, data = minDF),
  palette = "#2E9FDF",
  ggtheme = theme_minimal()
)

x <-
  anova(censored.clinical.res.cox,
        censored.clinical.biological.res.cox)
write.csv(
  as.data.frame(x),
  paste0(
    "data/survivalAnalysis/",
    clusterName,
    "DiseaseDurationClinicalVsBiologicalModel.csv"
  ),
  row.names = TRUE
)
print(x)
