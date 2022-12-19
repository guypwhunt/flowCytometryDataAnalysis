try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns
clusterNames <- clusterNames[4:3]
markersOrCells <- markersOrCellsClassification[2]

markersOrCell <- markersOrCells[1]

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    figureDirectory <- "data/pValueAdjustmentsResults/"

    fileNames <-
      list.files(
        figureDirectory,
        pattern = paste0("DifferentialStatesStatisticscsv", markersOrCell, ".csv")
      )

    fileNames <-
      fileNames[grepl(clusterName, fileNames, fixed = TRUE)]

    fileNamesPath <- paste0(figureDirectory, fileNames)

    dfs <- lapply(fileNamesPath, read.csv)

    names(dfs) <- fileNames

    i <- 1

    for (df in dfs) {
      df$experiment <- names(dfs)[i]
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

combinedDf <- combinedDf[combinedDf$fdr_adjusted_p_val < 0.05, ]

for (clusterName in clusterNames) {
  combinedDf <- combinedDf[combinedDf$typeOfCells %in% combinedDf[combinedDf$clusterModel == clusterName, "typeOfCells"],]
}

for(clusterName in clusterNames) {

  dsaResults <- combinedDf

  directoryNames <- unique(dsaResults$panel)
  directoryNames <- directoryNames[!is.na(directoryNames)]

  for (directory in directoryNames) {
    message(directory)

    dsResults <- dsaResults[dsaResults$panel == directory &
                              dsaResults$clusterModel == clusterName, ]

    clusterMedianValues <-
      read.csv(
        paste0(
          "data/",
          directory,
          "/clusteringOutput/",
          clusterName,
          markersOrCell,
          "Medians.csv"
        )
      )

    clusterMedianValues[clusterMedianValues$fileName ==  'BAS057_02', 'fileName'] <-
      'BAS_057_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00074-02', 'fileName'] <-
      'BLT00074-2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00075-04', 'fileName'] <-
      'BLT00075-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00092-06', 'fileName'] <-
      'BLT00092-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186-06', 'fileName'] <-
      'BLT00186-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186-09', 'fileName'] <-
      'BLT00186_9'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00211-03', 'fileName'] <-
      'BLT00211-3'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00211-06', 'fileName'] <-
      'BLT00211-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00211-06', 'fileName'] <-
      'BLT00211-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00214-02', 'fileName'] <-
      'BLT00214-2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00214-05', 'fileName'] <-
      'BLT00214-5'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00230-02', 'fileName'] <-
      'BLT00230-2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00243-07', 'fileName'] <-
      'BLT00243-7'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00244-04', 'fileName'] <-
      'BLT00244-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00244-06', 'fileName'] <-
      'BLT00244-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00244-06', 'fileName'] <-
      'BLT00244-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00254-05', 'fileName'] <-
      'BLT00254-5'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00254-07', 'fileName'] <-
      'BLT00254-7'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00265-04', 'fileName'] <-
      'BLT00265-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00265-04', 'fileName'] <-
      'BLT00265-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00265-08', 'fileName'] <-
      'BLT00265-8'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00271-04', 'fileName'] <-
      'BLT00271-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00271-07', 'fileName'] <-
      'BLT00271-7'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00274-06', 'fileName'] <-
      'BLT00274-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00285-03', 'fileName'] <-
      'BLT00285-3'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00285-05', 'fileName'] <-
      'BLT00285-5'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00285-05', 'fileName'] <-
      'BLT00285-5'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00286-04', 'fileName'] <-
      'BLT000286-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00286-06', 'fileName'] <-
      'BLT00286-6'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00297_02', 'fileName'] <-
      'BLT00297_2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT0244_02', 'fileName'] <-
      'BLT00244_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'QS_024-02', 'fileName'] <-
      'QS_024-2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00074-4', 'fileName'] <-
      'BLT00074_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00075-6', 'fileName'] <-
      'BLT00075_06'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186-2', 'fileName'] <-
      'BLT00186_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186-9', 'fileName'] <-
      'BLT00186_9'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00261-2', 'fileName'] <-
      'BLT00261_2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00297-2', 'fileName'] <-
      'BLT00297_2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00075-4_R1', 'fileName'] <-
      'BLT00075-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00211-3 ', 'fileName'] <-
      'BLT00211-3'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00297-2', 'fileName'] <-
      'BLT00297_2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT000214-5', 'fileName'] <-
      'BLT00214-5'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00057-4', 'fileName'] <-
      'BLT00057_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00092_4', 'fileName'] <-
      'BLT00092_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00198_2', 'fileName'] <-
      'BLT00198_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00214-4', 'fileName'] <-
      'BLT00214_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00230-4', 'fileName'] <-
      'BLT00230_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00242-2', 'fileName'] <-
      'BLT00242_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00243-5', 'fileName'] <-
      'BLT00243_05'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00274-2', 'fileName'] <-
      'BLT00274_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00243-5', 'fileName'] <-
      'BLT00243_05'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00274-2', 'fileName'] <-
      'BLT00274_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00274-4', 'fileName'] <-
      'BLT00274-05'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00286-5', 'fileName'] <-
      'BLT000286-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BAS057_02', 'fileName'] <-
      'BAS_057_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BAS033_02', 'fileName'] <-
      'BAS_033_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BAS057_02', 'fileName'] <-
      'BAS_057_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00074-4', 'fileName'] <-
      'BLT00074_04'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00075-6', 'fileName'] <-
      'BLT00075_06'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186_2', 'fileName'] <-
      'BLT00186_02'
    clusterMedianValues[clusterMedianValues$fileName ==  'BAS00101', 'fileName'] <-
      'BAS101'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00186-9', 'fileName'] <-
      'BLT00186_9'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00261_02', 'fileName'] <-
      'BLT00261_2'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00274-5', 'fileName'] <-
      'BLT00274-05'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00286-4', 'fileName'] <-
      'BLT000286-4'
    clusterMedianValues[clusterMedianValues$fileName ==  'BLT00297-2', 'fileName'] <-
      'BLT00297_2'

    clusterMedianValues <-
      clusterMedianValues[clusterMedianValues[, clusterName] %in% unique(dsResults$cluster_id), ]

    combinedDf <-
      merge(
        clusterMedianValues[, c("fileName", clusterName, "GPR32_median")],
        dsResults[, c("cluster_id", "typeOfCells")],
        by.x = clusterName,
        by.y = "cluster_id",
        all.x = TRUE
      )

    combinedDf$cell_population_cluster <-
      #paste(combinedDf$typeOfCells, combinedDf[, clusterName], sep = " ")
      combinedDf$typeOfCells

    combinedDf$cell_population_cluster <-
      str_replace_all(combinedDf$cell_population_cluster, " ", "_")

    combinedDf <- distinct(combinedDf)

    combinedDfT <-
      reshape(
        combinedDf[, c("fileName", "cell_population_cluster", "GPR32_median")],
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
        directory,
        clusterName,
        markersOrCell,
        "_medianValues.csv"
      ),
      row.names = FALSE
    )
  }
}
