try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

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

combinedDf <- combinedDf[!is.na(combinedDf$outcomeDeathDate), ]

combinedDf$fileName <- factor(combinedDf$fileName)

combinedDf <- combinedDf[combinedDf$visit == 1, ]

scaledDF <- combinedDf
columnNames <- c(
  colnames(clusterMedianValues)[
    colnames(clusterMedianValues) != "fileName"],
  "survivalFromVisit",
  "alsfrsR",
  "timeFromOnsetToVisitInYears",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "ageAtVisit"
)

scaledDF[, columnNames] <- scale(scaledDF[, columnNames])

columnNames <- columnNames[columnNames != "survivalFromVisit"]

for (column in columnNames) {
  print(plot(density(na.omit(scaledDF[, column])), main = column))
  print(hist(na.omit(scaledDF[, column]), main = column))
  print(column)
  print(cor(scaledDF$survivalFromVisit, scaledDF[, column]))
  plot(as.formula(paste0("survivalFromVisit ~ ", column)), data = scaledDF,
       main = column)
  }

# https://stackoverflow.com/questions/33740012/command-for-finding-the-best-linear-model-in-r
library(MASS)
fitx <- lm(
  as.formula(paste("survivalFromVisit ~ 0",
                   sapply(list(columnNames), paste, collapse=" + "),

                   sep=" + "))
  ,
  data = scaledDF
)
step <- stepAIC(fitx, direction = "both")
step$anova # display results

fit <-
  lm(
    survivalFromVisit ~ GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2 +
      GPR32_median.monocytes_Clusters_meta_clusters_flowsom_4 +
      timeFromOnsetToVisitInYears + diagnosticDelayInYears + ageAtOnset -
      1
    ,
    data = scaledDF
  )

summary(fit)
plot(fit)

library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
pred.mm <-
  ggpredict(fit, terms = c("GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2"))  # this gives overall predictions for the model
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
      data = scaledDF,
      # adding the raw data (scaled values)
      aes(x = GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2, y = survivalFromVisit, colour = fastSlow)
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
      data = scaledDF,
      # adding the raw data (scaled values)
      aes(x = GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2, y = survivalFromVisit, colour = ethnicity)
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
      data = scaledDF,
      # adding the raw data (scaled values)
      aes(x = GPR32_median.monocytes_Clusters_meta_clusters_flowsom_2, y = survivalFromVisit, colour = BulbarLimb)
    ) +
    labs(x = "Median GPR32 Abundance", y = "Survival from Visit in Years") +
    guides(color = guide_legend(title = "Site of Onset")) +
    scale_colour_viridis_d() +
    theme_minimal()
)

(
  ggplot(data = scaledDF,
         # adding the raw data (scaled values)
         aes(x = ageAtVisit, y = survivalFromVisit)) +
    geom_point() +
    labs(x = "Age at Visit", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

(
  ggplot(scaledDF, aes(x = alsfrsR, y = survivalFromVisit)) +
    geom_point() +
    labs(x = "ALSFRSR", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

(
  ggplot(
    scaledDF,
    aes(x = timeFromOnsetToVisitInYears, y = survivalFromVisit)
  ) +
    geom_point() +
    labs(x = "Time from Onset to Visit", y = "Survival from Visit in Years") +
    geom_smooth(method = "lm") +
    theme_minimal()
)

boxplot(scaledDF$survivalFromVisit ~ scaledDF$gender)
boxplot(scaledDF$survivalFromVisit ~ scaledDF$ethnicity)
boxplot(scaledDF$survivalFromVisit ~ scaledDF$BulbarLimb)
