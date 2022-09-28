try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library("survival")
library("survminer")

loadlibraries()

results <- read.csv("data/medianValues/medianValues.csv")

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo[which(experimentInfo[, "sample_id"] == "BAS_057_02", arr.ind =
                       TRUE), "sample_id"] <- "BAS057_02"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00074-2", arr.ind =
                       TRUE), "sample_id"] <- "BLT00074-02"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00075-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00075-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00092-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00092-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00186-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_9", arr.ind =
                       TRUE), "sample_id"] <- "BLT00186-09"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-3", arr.ind =
                       TRUE), "sample_id"] <- "BLT00211-03"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00211-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00211-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214-2", arr.ind =
                       TRUE), "sample_id"] <- "BLT00214-02"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214-5", arr.ind =
                       TRUE), "sample_id"] <- "BLT00214-05"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00230-2", arr.ind =
                       TRUE), "sample_id"] <- "BLT00230-02"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00243-7", arr.ind =
                       TRUE), "sample_id"] <- "BLT00243-07"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00244-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00244-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00244-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00254-5", arr.ind =
                       TRUE), "sample_id"] <- "BLT00254-05"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00254-7", arr.ind =
                       TRUE), "sample_id"] <- "BLT00254-07"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00265-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00265-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-8", arr.ind =
                       TRUE), "sample_id"] <- "BLT00265-08"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00271-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00271-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00271-7", arr.ind =
                       TRUE), "sample_id"] <- "BLT00271-07"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00274-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-3", arr.ind =
                       TRUE), "sample_id"] <- "BLT00285-03"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-5", arr.ind =
                       TRUE), "sample_id"] <- "BLT00285-05"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-5", arr.ind =
                       TRUE), "sample_id"] <- "BLT00285-05"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT000286-4", arr.ind =
                       TRUE), "sample_id"] <- "BLT00286-04"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00286-6", arr.ind =
                       TRUE), "sample_id"] <- "BLT00286-06"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00297_2", arr.ind =
                       TRUE), "sample_id"] <- "BLT00297_02"
experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244_02", arr.ind =
                       TRUE), "sample_id"] <- "BLT0244_02"
experimentInfo[which(experimentInfo[, "sample_id"] == "QS_024-2", arr.ind =
                       TRUE), "sample_id"] <- "QS_024-02"

combinedDf <-
  merge(
    results,
    experimentInfo,
    by.x = "fileName",
    by.y = "sample_id",
    all.x = TRUE
  )

#combinedDf <- combinedDf[!is.na(combinedDf$outcomeDeathDate), ]

dateOfLateOutCome <- combinedDf[is.na(combinedDf$outcomeDeathDate), c("outcomeLastVisitDate", "visitDate")]

dateOfLateOutCome <- (dateOfLateOutCome$outcomeLastVisitDate - dateOfLateOutCome$visitDate)/(24*60*60*365.25)

combinedDf[is.na(combinedDf$outcomeDeathDate), "survivalFromVisit"] <- dateOfLateOutCome

combinedDf <- combinedDf[combinedDf$caseControl == "Case", ]

combinedDf$fileName <- factor(combinedDf$fileName)

combinedDf <- combinedDf[combinedDf$visit == 1, ]

combinedDf$survivalFromVisit <- combinedDf$survivalFromVisit*365.25

combinedDf$outcome <- factor(combinedDf$outcome)
combinedDf$status <- as.numeric(combinedDf$outcome)

combinedDf$gender <- factor(combinedDf$gender)
combinedDf$sex <- as.numeric(combinedDf$gender)

combinedDf$ethnicity <- factor(combinedDf$ethnicity)
combinedDf$ethnicityID <- as.numeric(combinedDf$ethnicity)

combinedDf$BulbarLimb <- factor(combinedDf$BulbarLimb)
combinedDf$onset <- as.numeric(combinedDf$BulbarLimb)


res.cox <- coxph(Surv(survivalFromVisit, status) ~ gender, data = combinedDf)


covariates <- c("medianValue",
                  "ageAtVisit",
                  "sex",
                  "ethnicityID",
                  "onset",
                  "alsfrsR",
                  "timeFromOnsetToVisitInYears",
                  "diagnosticDelayInYears",
                  "ageAtOnset")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(survivalFromVisit, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = combinedDf)})
# Extract data
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

res.cox <- coxph(Surv(survivalFromVisit, status) ~ ageAtVisit + ageAtOnset + medianValue, data = combinedDf)
summary(res.cox)

ggsurvplot(survfit(res.cox, data = combinedDf), palette = "#2E9FDF",
           ggtheme = theme_minimal())

sex_df <- with(combinedDf,
               data.frame(sex = c(1, 2),
                          medianValue = rep(mean(medianValue, na.rm = TRUE), 2),
                          ageAtOnset = rep(mean(ageAtOnset, na.rm = TRUE), 2),
                          ageAtVisit = rep(mean(ageAtVisit, na.rm = TRUE), 2),
                          diagnosticDelayInYears = rep(mean(diagnosticDelayInYears, na.rm = TRUE), 2)
               )
)
sex_df

fit <- survfit(res.cox, data = combinedDf, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())


library(StepReg)
#install.packages("StepReg")

stepwiseCox(
  Surv(survivalFromVisit) ~ ageAtVisit + ageAtOnset + medianValue,
  data = combinedDf,
  selection = "bidirection"
)

library(My.stepwise)
#install.packages("My.stepwise")

minDF <- combinedDf[, c("survivalFromVisit", "status", covariates)]
minDF <- na.omit(minDF)
minDF$status <- as.double(minDF$status)
minDF$sex <- as.factor(minDF$sex)
minDF$ethnicityID <- as.factor(minDF$ethnicityID)
minDF$onset <- as.factor(minDF$onset)

My.stepwise.coxph(Time = "survivalFromVisit", Status = "status",
                  variable.list = covariates[c(1,3)], data = minDF)

for (covariate in covariates){
  print(covariate)
  print(typeof(minDF[,covariate]))
}
