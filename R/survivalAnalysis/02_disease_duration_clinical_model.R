try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo <- experimentInfo[experimentInfo$visit == 1,]
experimentInfo <-
  experimentInfo[experimentInfo$experiment == "flowCytometry",]
experimentInfo <-
  experimentInfo[experimentInfo$caseControl == "Case",]

combinedDf <-
  experimentInfo

dateOfLateOutCome <-
  combinedDf[is.na(combinedDf$outcomeDeathDate), c("outcomeLastVisitDate", "visitDate")]

dateOfLateOutCome <-
  (dateOfLateOutCome$outcomeLastVisitDate - dateOfLateOutCome$visitDate) /
  (24 * 60 * 60 * 365.25)

combinedDf[is.na(combinedDf$outcomeDeathDate), "diseaseDurationInYears"] <-
  dateOfLateOutCome

combinedDf$fileName <- factor(combinedDf$patient_id)

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
  "sex",
  "ethnicityID",
  "onset",
  "alsfrsR",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "riluzole",
  "delataAlsfrsRScore"
)

combinedDf$sex <- as.factor(combinedDf$sex)
combinedDf$ethnicityID <- as.factor(combinedDf$ethnicityID)
combinedDf$onset <- as.factor(combinedDf$onset)

# Censored Analysis
minDF <-
  combinedDf[, c("diseaseDurationInYears",
                 "status",
                 clinicalCovariates)]
minDF <- na.omit(minDF)

# Univariate Analysis
univ_formulas <- sapply(clinicalCovariates,
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
res <- as.data.frame(res)

write.csv(res,
          "data/survivalAnalysis/clinicalUnivariateAnalysis.csv")

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
    Surv(diseaseDurationInYears, status) ~ ethnicityID + alsfrsR +
      ageAtOnset + riluzole + delataAlsfrsRScore - 1,
    data = minDF
  )
censored.clinical.res.cox.summary <-
  summary(censored.clinical.res.cox)
saveRDS(censored.clinical.res.cox, file = paste0(
  "data/survivalAnalysis/",
  "DiseaseDurationClinicalModel.rds"
))

x <- as.data.frame(censored.clinical.res.cox.summary$coefficients)

colnames(x) <- c( "Regression Coefficient", "Hazard Ratio",
                 "Confidence Interval", "Statistical Significance",
                 "Global Significance (P-Value)")


write.csv(
  x,
  paste0(
    "data/survivalAnalysis/",
    "DiseaseDurationClinicalModel.csv"
  ),
  row.names = TRUE
)

ggsurvplot(
  survfit(censored.clinical.res.cox, data = minDF),
  palette = "#2E9FDF",
  ggtheme = theme_minimal()
)

rm(list=ls())

