try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

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

combinedDf <- combinedDf[!is.na(combinedDf$outcomeDeathDate),]

combinedDf$fileName <- factor(combinedDf$fileName)

combinedDf <- combinedDf[combinedDf$visit == 1,]

minMaxDf <- combinedDf
columnNames <- c(
  "medianValue",
  "survivalFromVisit",
  "alsfrsR",
  "timeFromOnsetToVisitInYears",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "ageAtVisit"
)

minMaxDf[, columnNames] <- scale(minMaxDf[, columnNames])

print(plot(density(minMaxDf[, columnNames[1]])))
print(plot(density(minMaxDf[, columnNames[2]])))

#for (col in columnNames) {
#  dfMax <- max(minMaxDf[, col])
#  dfMin <- min(minMaxDf[, col])

#  scaledColumn <- (minMaxDf[, col] - dfMin) / (dfMax - dfMin)

#  minMaxDf[, col] <- scaledColumn

#  d <- density(scaledColumn)
#  print(plot(d))
#}


# https://www.scribbr.com/statistics/linear-regression-in-r/
# Initial Checks
# In the data normally distributed?
hist(minMaxDf$medianValue)
hist(minMaxDf$survivalFromVisit)

# Are the relationships linear?
cor(minMaxDf$survivalFromVisit, minMaxDf$medianValue)
plot(survivalFromVisit ~ medianValue, data = minMaxDf)

# https://stackoverflow.com/questions/33740012/command-for-finding-the-best-linear-model-in-r
library(MASS)
fitx <- lm(
  survivalFromVisit ~ 0 +
    medianValue +
    ageAtVisit +
    gender +
    ethnicity +
    BulbarLimb +
    alsfrsR +
    timeFromOnsetToVisitInYears +
    diagnosticDelayInYears +
    ageAtOnset
  ,
  data = minMaxDf
)
step <- stepAIC(fitx, direction = "both")
step$anova # display results

fit <-
  lm(
    survivalFromVisit ~ medianValue + ageAtVisit + gender + ethnicity +
      BulbarLimb + alsfrsR + timeFromOnsetToVisitInYears - 1
    ,
    data = minMaxDf
  )

summary(fit)
plot(fit)

#boxplot(minMaxDf$survivalFromVisit ~ minMaxDf$gender)


library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
pred.mm <-
  ggpredict(fit, terms = c("medianValue"))  # this gives overall predictions for the model
# Plot the predictions

(
  ggplot(pred.mm) +
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(
      aes(
        x = x,
        ymin = predicted - std.error,
        ymax = predicted + std.error
      ),
      fill = "lightgrey",
      alpha = 0.5
    ) +  # error band
    geom_point(
      data = minMaxDf,
      # adding the raw data (scaled values)
      aes(x = medianValue, y = survivalFromVisit, colour = fastSlow)
    ) +
    labs(x = "Median GPR32 Abundance", y = "Survival from Visit in Years") +
    guides(color = guide_legend(title = "Progression")) +
    scale_colour_viridis_d() +
    theme_minimal()
)

(
  ggplot(pred.mm) +
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(
      aes(
        x = x,
        ymin = predicted - std.error,
        ymax = predicted + std.error
      ),
      fill = "lightgrey",
      alpha = 0.5
    ) +  # error band
    geom_point(
      data = minMaxDf,
      # adding the raw data (scaled values)
      aes(x = medianValue, y = survivalFromVisit, colour = ethnicity)
    ) +
    labs(x = "Median GPR32 Abundance", y = "Survival from Visit in Years") +
    guides(color = guide_legend(title = "Ethnicity")) +
    scale_colour_viridis_d() +
    theme_minimal()
)

(
  ggplot(pred.mm) +
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(
      aes(
        x = x,
        ymin = predicted - std.error,
        ymax = predicted + std.error
      ),
      fill = "lightgrey",
      alpha = 0.5
    ) +  # error band
    geom_point(
      data = minMaxDf,
      # adding the raw data (scaled values)
      aes(x = medianValue, y = survivalFromVisit, colour = BulbarLimb)
    ) +
    labs(x = "Median GPR32 Abundance", y = "Survival from Visit in Years") +
    guides(color = guide_legend(title = "Site of Onset")) +
    scale_colour_viridis_d() +
    theme_minimal()
)

(
  ggplot(data = minMaxDf,
         # adding the raw data (scaled values)
         aes(x = ageAtVisit, y = survivalFromVisit)) +
    geom_point() +
    labs(x = "Age at Visit", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

(
  ggplot(minMaxDf, aes(x = alsfrsR, y = survivalFromVisit)) +
    geom_point() +
    labs(x = "ALSFRSR", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

(
  ggplot(
    minMaxDf,
    aes(x = timeFromOnsetToVisitInYears, y = survivalFromVisit)
  ) +
    geom_point() +
    labs(x = "Time from Onset to Visit", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

boxplot(minMaxDf$survivalFromVisit ~ minMaxDf$gender)
boxplot(minMaxDf$survivalFromVisit ~ minMaxDf$ethnicity)
boxplot(minMaxDf$survivalFromVisit ~ minMaxDf$BulbarLimb)


anova(mixedFit)


# Linear mixed model
# https://ourcodingclub.github.io/tutorials/mixed-models/
library(lme4)
mixedFit <-
  lmer(
    survivalFromVisit ~ medianValue + ageAtVisit + (1 |
                                                      ethnicity) + (1 |
                                                                      BulbarLimb) +
      alsfrsR + timeFromOnsetToVisitInYears - 1,
    data = combinedDf
  )
summary(mixedFit)
coef(mixedFit)
plot(mixedFit)
qqnorm(resid(mixedFit))
qqline(resid(mixedFit))
