try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    medianResultsAllClusters <- read.csv("data/medianValues/medianValues.csv")

    medianResultsAllClusters <- medianResultsAllClusters[
      medianResultsAllClusters$clusterName == clusterName &
      medianResultsAllClusters$markersOrCell == markersOrCell,]

    head(medianResultsAllClusters)

    medianResultsIndividualClusters <-
      data.frame(fileName = as.character(),
                 clusterName = as.character(),
                 clusterId = as.character(),
                 markersOrCell = as.character(),
                 GPR32_median = as.double()
                 )

    for (panel in unique(medianResultsAllClusters$panel)){
      x <-
        read.csv(paste0("data/monocytes/clusteringOutput/",clusterName, markersOrCell, "Medians.csv"))

      x$clusterName <- clusterName
      x$markersOrCell <- markersOrCell
      x$clusterId <- x[, clusterName]

      medianResultsIndividualClusters <- rbind(medianResultsIndividualClusters, x)

    }
  }
}


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
    medianResultsAllClusters,
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

unique(dsResults$panel)

head(monocytesResults)
head(dsResults)

monocytesDsResults <- dsResults[dsResults$panel == "monocytes"]

#monocytesResults$meta_clusters_flowsom  monocytesDsResults$cluster_id

monocytesResults <-
  monocytesResults[monocytesResults[, 'meta_clusters_flowsom'] %in% monocytesDsResults$cluster_id,]

monocytesResults <-
  monocytesResults[, c("fileName", "meta_clusters_flowsom", "GPR32_median")]

monocytesResults <-
  reshape(
    monocytesResults,
    idvar = "fileName",
    timevar = "meta_clusters_flowsom",
    direction = "wide"
  )

combinedDf <- merge(combinedDf,
                    as.data.frame(monocytesResults),
                    by = "fileName",
                    all.x = TRUE
                    )

colnames(combinedDf)[colnames(combinedDf)== "medianValue"] <- "GPR32_median.all"

combinedDf <- combinedDf[!is.na(combinedDf$outcomeDeathDate), ]

combinedDf$fileName <- factor(combinedDf$fileName)

combinedDf <- combinedDf[combinedDf$visit == 1, ]

scaledDF <- combinedDf
columnNames <- c(
  "GPR32_median.all",
  "GPR32_median.11",
  "GPR32_median.3",
  "GPR32_median.2",
  "GPR32_median.4",
  "survivalFromVisit",
  "alsfrsR",
  "timeFromOnsetToVisitInYears",
  "diagnosticDelayInYears",
  "ageAtOnset",
  "ageAtVisit"
)

scaledDF[, columnNames] <- scale(scaledDF[, columnNames])

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
                   sapply(list(columnNames[columnNames!= "survivalFromVisit"]), paste, collapse=" + "),

                   sep=" + "))
  ,
  data = scaledDF
)
step <- stepAIC(fitx, direction = "both")
step$anova # display results

fit <-
  lm(
    survivalFromVisit ~ GPR32_median.2 + GPR32_median.4 + timeFromOnsetToVisitInYears +
      diagnosticDelayInYears + ageAtOnset - 1
    ,
    data = scaledDF
  )

summary(fit)
plot(fit)

library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
pred.mm <-
  ggpredict(fit, terms = c("GPR32_median.2"))  # this gives overall predictions for the model
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
      aes(x = GPR32_median.2, y = survivalFromVisit, colour = fastSlow)
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
      aes(x = GPR32_median.2, y = survivalFromVisit, colour = ethnicity)
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
      aes(x = GPR32_median.2, y = survivalFromVisit, colour = BulbarLimb)
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
