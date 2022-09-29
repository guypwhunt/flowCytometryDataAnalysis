try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))
library(survival)
library(survminer)

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification
directoryNames <- c("bCells", "monocytes", "tCells", "senescence")

clusterName <- clusterNames[3]
markersOrCell <- markersOrCells[1]
panel <- directory <- directoryNames[2]

clusterMedianValues <- read.csv(paste0("data/medianValues/",directory, clusterName, markersOrCell, "_medianValues.csv"))

#results <- read.csv("data/medianValues/medianValues.csv")

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

combinedDf <-
  merge(
    clusterMedianValues,
    experimentInfo,
    by.x = "fileName",
    by.y = "sample_id",
    all.x = TRUE
  )

dsResults <- read.csv(
  paste0(
    "data/pValueAdjustmentsResults/",
    clusterName,
    "fastSlowVisits1AllCellsDifferentialStatesStatisticscsv",
    markersOrCell ,
    ".csv"
  )
)

dsResults <- dsResults[dsResults$fdr_adjusted_p_val <= 0.05,]

dsResults <- dsResults[dsResults$panel == directory]

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
combinedDf$status <- ifelse(combinedDf$status==2,1,0)

combinedDf$gender <- factor(combinedDf$gender)
combinedDf$sex <- as.numeric(combinedDf$gender)

combinedDf$ethnicity <- factor(combinedDf$ethnicity)
combinedDf$ethnicityID <- as.numeric(combinedDf$ethnicity)

combinedDf$BulbarLimb <- factor(combinedDf$BulbarLimb)
combinedDf$onset <- as.numeric(combinedDf$BulbarLimb)

covariates <- c(  colnames(clusterMedianValues)[
  colnames(clusterMedianValues) != "fileName"],
                "ageAtVisit",
                "sex",
                "ethnicityID",
                "onset",
                "alsfrsR",
                "timeFromOnsetToVisitInYears",
                "diagnosticDelayInYears",
                "ageAtOnset")

combinedDf$sex <- as.factor(combinedDf$sex)
combinedDf$ethnicityID <- as.factor(combinedDf$ethnicityID)
combinedDf$onset <- as.factor(combinedDf$onset)

deadMinDF <- combinedDf
deadMinDF <- deadMinDF[deadMinDF$outcome == "dead", ]
deadMinDF <- deadMinDF[, c("survivalFromVisit", "status", covariates)]
deadMinDF <- na.omit(deadMinDF)
deadMinDF$ethnicityID <- factor(deadMinDF$ethnicityID)

# Univariate Analysis
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(survivalFromVisit, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = deadMinDF)})
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

univ_results$ethnicityID <- univ_results$ethnicityID[c(1,2,4,5)]
names(univ_results$ethnicityID) <- names(univ_results$sex)
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

colnames(deadMinDF)[
  !colnames(deadMinDF) %in% c("survivalFromVisit", "status")]

res.cox <-
  coxph(as.formula(paste(
    "Surv(survivalFromVisit, status) ~ 0",
    sapply(list(covariates), paste, collapse =
             " + "),

    sep = " + "
  ))
  , data = deadMinDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

res.cox <- coxph(Surv(survivalFromVisit, status) ~ GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2 +
                   GPR32_median.monocytes_Clusters_meta_clusters_flowsom_4 +
                   ageAtVisit + onset + alsfrsR - 1, data = deadMinDF)
summary(res.cox)

ggsurvplot(survfit(res.cox, data = deadMinDF), palette = "#2E9FDF",
           ggtheme = theme_minimal())


# Censored Analysis
minDF <- combinedDf[, c("survivalFromVisit", "status", covariates)]
minDF <- na.omit(minDF)

# Univariate Analysis
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(survivalFromVisit, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = minDF)})
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

univ_results$ethnicityID <- univ_results$ethnicityID[c(1,2,4,5)]
names(univ_results$ethnicityID) <- names(univ_results$sex)
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

res.cox <-
  coxph(as.formula(paste(
    "Surv(survivalFromVisit, status) ~ 0",
    sapply(list(covariates), paste, collapse =
             " + "),

    sep = " + "
  ))
  , data = minDF)

library(MASS)
step <- stepAIC(res.cox, direction = "both")
step$anova

res.cox <- coxph(Surv(survivalFromVisit, status) ~ ageAtVisit + onset + alsfrsR -
                   1, data = minDF)
summary(res.cox)

ggsurvplot(survfit(res.cox, data = minDF), palette = "#2E9FDF",
           ggtheme = theme_minimal())
