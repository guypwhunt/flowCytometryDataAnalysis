try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo <- experimentInfo[experimentInfo$visit == 1, ]
#experimentInfo <-
#  experimentInfo[experimentInfo$experiment == "lipidomics", ]
experimentInfo <-
  experimentInfo[experimentInfo$caseControl == "Case", ]

mergedDF <- read.csv(
  "data/lipidomics/normalisedasinhTransformedExpressionDataRawOutliersAndDuplicatesRemoved.csv",
  row.names = 1
)

mergedDF <- mergedDF[rownames(mergedDF) %in% c("RvD1", "PGD2", "15R-LXB4", "7S,14S-diHDHA"), ]

mergedDF <- t(mergedDF)
mergedDF <- as.data.frame(mergedDF)
mergedDF$classification <- row.names(mergedDF)

colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "-", "")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), " ", "")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), ",", "")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "1", "one")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "2", "two")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "3", "three")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "4", "four")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "5", "five")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "6", "six")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "7", "seven")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "8", "eight")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "9", "nine")
colnames(mergedDF) <- str_replace_all(colnames(mergedDF), "0", "zero")

combinedDf <- experimentInfo

#combinedDf <- merge(
#  experimentInfo,
#  mergedDF,
#  by.x = "classification",
#  by.y = "classification",
#  all.x = TRUE
#)

dateOfLateOutCome <-
  combinedDf[is.na(combinedDf$outcomeDeathDate), c("outcomeLastVisitDate", "visitDate")]

dateOfLateOutCome <-
  (dateOfLateOutCome$outcomeLastVisitDate - dateOfLateOutCome$visitDate) /
  (24 * 60 * 60 * 365.25)

combinedDf[is.na(combinedDf$outcomeDeathDate), "diseaseDurationInYears"] <-
  dateOfLateOutCome

combinedDf$diseaseDurationInYears <-
  combinedDf$diseaseDurationInYears* 365.25

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
  "sex",
  "ethnicityID",
  "onset",
  "alsfrsR",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "delataAlsfrsRScore"
)

biologicalCovariates <-
  c(colnames(mergedDF)[colnames(mergedDF) != "classification"])

clinicalAndBiologicalCovariates <-
  append(clinicalCovariates, biologicalCovariates)

combinedDf$sex <- as.factor(combinedDf$sex)
combinedDf$ethnicityID <- as.factor(combinedDf$ethnicityID)
combinedDf$onset <- as.factor(combinedDf$onset)

# Censored Analysis
minDF <-
  #combinedDf[, c("diseaseDurationInYears",
  #               "status",
  #               clinicalAndBiologicalCovariates)]
  combinedDf[, c("diseaseDurationInYears",
                 "status",
                 clinicalCovariates)]

minDF <- minDF[rowSums(is.na(minDF)) < 1,]

clinicalAndBiologicalCovariates <- colnames(minDF)[!colnames(minDF) %in% c("diseaseDurationInYears",
                                                                          "status")]


# Univariate Analysis
univ_formulas <- sapply(clinicalAndBiologicalCovariates,
                        function(x)
                          as.formula(paste('Surv(diseaseDurationInYears, status)~', x)))

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

univ_results$ethnicityID <-
  univ_results$ethnicityID[c(1, 2, 4, 5)]
names(univ_results$ethnicityID) <- names(univ_results$sex)
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

res.cox <-
  coxph(as.formula(paste(
    "Surv(diseaseDurationInYears, status) ~ 0",
    sapply(list(clinicalCovariates), paste, collapse =
             " + "),

    sep = " + "
  ))
  , data = minDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

censored.clinical.res.cox <-
  coxph(
    Surv(diseaseDurationInYears, status) ~ 0 + sex + ethnicityID +
      onset + alsfrsR + diagnosticDelayInYears + ageAtOnset + delataAlsfrsRScore,
    data = minDF
  )
saveRDS(censored.clinical.res.cox, file = paste0(
  "data/survivalAnalysis/",
  "lipidomics",
  "diseaseDurationInYearsClinicalModel.rds"
))

censored.clinical.res.cox.summary <-
  summary(censored.clinical.res.cox)
write.csv(
  as.data.frame(censored.clinical.res.cox.summary$coefficients),
  paste0(
    "data/survivalAnalysis/",
    "lipidomics",
    "diseaseDurationInYearsClinicalModel.csv"
  ),
  row.names = TRUE
)

ggsurvplot(
  survfit(censored.clinical.res.cox, data = minDF),
  palette = "#2E9FDF",
  ggtheme = theme_minimal()
)

clinicalCovariates <- clinicalCovariates[clinicalCovariates %in% c("sex",  "ethnicityID",
                                                                     "onset", "alsfrsR", "diagnosticDelayInYears", "ageAtOnset", "delataAlsfrsRScore")]
biologicalCovariates <- biologicalCovariates[biologicalCovariates %in% colnames(mergedDF)]

clinicalAndBiologicalCovariates <-
  append(clinicalCovariates, biologicalCovariates)

res.cox <-
  coxph(as.formula(paste(
    "Surv(diseaseDurationInYears, status) ~ 0",
    sapply(list(clinicalAndBiologicalCovariates), paste, collapse =
             " + "),

    sep = " + "
  ))
  , data = minDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

censored.clinical.biological.res.cox <-
  coxph(
    Surv(diseaseDurationInYears, status) ~ sex + ethnicityID + onset +
      alsfrsR + diagnosticDelayInYears + ageAtOnset + delataAlsfrsRScore +
      RvDone + sevenSonefourSdiHDHA - 1,
    data = minDF
  )

saveRDS(censored.clinical.biological.res.cox, file = paste0(
  "data/survivalAnalysis/",
  "lipidomics",
  "diseaseDurationInYearsBiologicalModel.rds"
))
censored.clinical.biological.res.cox.summary <-
  summary(censored.clinical.biological.res.cox)
write.csv(
  as.data.frame(censored.clinical.biological.res.cox.summary$coefficients),
  paste0(
    "data/survivalAnalysis/",
    "lipidomics",
    "diseaseDurationInYearsBiologicalModel.csv"
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
    "lipidomics",
    "diseaseDurationInYearsClinicalVsBiologicalModel.csv"
  ),
  row.names = TRUE
)
print(x)
