try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns
clusterNames <- clusterNames[4:3]
markersOrCells <- markersOrCellsClassification[3]
markerNames <- c("GPR18", "GPR32")

markersOrCell <- markersOrCells[1]
clusterName <- clusterNames[1]
markerName <- markerNames[1]

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

for (markerName in markerNames) {
  for (clusterName in clusterNames) {
    for (markersOrCell in markersOrCells) {
      figureDirectory <- "data/pValueAdjustmentsResults/"

      fileNames <-
        list.files(
          figureDirectory,
          pattern = paste0(
            clusterName,
            markersOrCell,
            "DifferentialStatesStatisticscsv.csv"
          )
        )

      fileNamesPath <- paste0(figureDirectory, fileNames)

      dfs <- lapply(fileNamesPath, read.csv)

      names(dfs) <- fileNames

      i <- 1

      for (df in dfs) {
        df$clusterModel <- clusterName
        i <- i + 1

        if (exists("combinedDf")) {
          combinedDf <- rbind(combinedDf, df)
        } else {
          combinedDf <- df
        }
      }
    }
  }

  combinedDf <- combinedDf[combinedDf$Marker == markerName,]
  combinedDf <- combinedDf[combinedDf$FDR.Adjusted.P.Value < 0.05,]

  for (clusterName in clusterNames) {
    combinedDf <-
      combinedDf[combinedDf$Cell.Population.Name %in% combinedDf[combinedDf$clusterModel == clusterName, "Cell.Population.Name"],]
  }

  for (clusterName in clusterNames) {
    message(clusterName)

    dsaResults <- combinedDf

    directoryNames <- unique(dsaResults$Panel)
    directoryNames <- directoryNames[!is.na(directoryNames)]

    for (directory in directoryNames) {
      message(directory)

      dsResults <- dsaResults[dsaResults$Panel == directory &
                                dsaResults$clusterModel == clusterName,]

      clusterMedianValues <-
        read.csv(
          paste0(
            "data/",
            directory,
            "/clusteringOutput/",
            markerName,
            clusterName,
            markersOrCell,
            "Medians.csv"
          )
        )

      clusterMedianValues <- updateDfFileNames(clusterMedianValues)

      clusterMedianValues <-
        clusterMedianValues[clusterMedianValues[, clusterName] %in% unique(dsResults$ID),]

      combinedDf <-
        merge(
          clusterMedianValues[, c("fileName", clusterName, paste0(markerName, "_median"))],
          dsResults[, c("ID", "Cell.Population.Name")],
          by.x = clusterName,
          by.y = "ID",
          all.x = TRUE
        )

      combinedDf$cell_population_cluster <-
        combinedDf$Cell.Population.Name

      combinedDf$cell_population_cluster <-
        str_replace_all(combinedDf$cell_population_cluster, " ", "_")

      combinedDf <- distinct(combinedDf)

      combinedDfT <-
        reshape(
          combinedDf[, c("fileName", "cell_population_cluster", paste0(markerName, "_median"))],
          idvar = "fileName",
          timevar = "cell_population_cluster",
          direction = "wide"
        )

      colnames(combinedDfT) <-
        str_replace_all(colnames(combinedDfT), "/", "_Or_")

      colnames(combinedDfT) <-
        str_replace_all(colnames(combinedDfT), "\\+", "_Positive_")

      colnames(combinedDfT) <-
        str_replace_all(colnames(combinedDfT), "-", "_Negative_")

      colnames(combinedDfT) <-
        str_replace_all(colnames(combinedDfT), "\\.", "_Positive_")

      write.csv(
        combinedDfT,
        paste0(
          "data/medianValues/",
          markerName,
          directory,
          clusterName,
          markersOrCell,
          "_medianValues.csv"
        ),
        row.names = FALSE
      )
    }
  }


  rm(list = ls()[!ls() %in% c("clusterNames", "markersOrCells", "markerNames",
                                "clusterName", "markersOrCell", "markerName",
                                "experimentInfo")])
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))
}
