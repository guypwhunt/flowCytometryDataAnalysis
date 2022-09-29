try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification
directoryNames <- c("bCells", "monocytes", "tCells", "senescence")

clusterName <- clusterNames[3]
markersOrCell <- markersOrCells[1]
panel <- directory <- directoryNames[2]

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)
for (directory in directoryNames) {
  for (clusterName in clusterNames) {
    for (markersOrCell in markersOrCells) {
      dsResults <- read.csv(
        paste0(
          "data/pValueAdjustmentsResults/",
          clusterName,
          "fastSlowVisits1AllCellsDifferentialStatesStatisticscsv",
          markersOrCell ,
          ".csv"
        )
      )

      dsResults <- dsResults[dsResults$fdr_adjusted_p_val <= 0.05, ]

      dsResults <- dsResults[dsResults$panel == directory, ]

      medianResultsAllClusters <-
        read.csv("data/medianValues/medianValues.csv")

      medianResultsAllClusters <- medianResultsAllClusters[medianResultsAllClusters$clusterName == clusterName &
                                                             medianResultsAllClusters$markersOrCell == markersOrCell &
                                                             medianResultsAllClusters$panel == directory, ]

      head(medianResultsAllClusters)

      medianResultsIndividualClusters <-
        data.frame(
          fileName = as.character(),
          panel = as.character(),
          clusterName = as.character(),
          clusterId = as.character(),
          markersOrCell = as.character(),
          GPR32_median = as.double()
        )

        x <-
          read.csv(
            paste0(
              "data/",
              panel,
              "/clusteringOutput/",
              clusterName,
              markersOrCell,
              "Medians.csv"
            )
          )

        x$panel <- panel
        x$clusterName <- clusterName
        x$markersOrCell <- markersOrCell
        x$clusterId <- x[, clusterName]

        x <- x[, colnames(medianResultsIndividualClusters)]

        medianResultsIndividualClusters <-
          rbind(medianResultsIndividualClusters, x)

      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BAS057_02', 'fileName'] <-
        'BAS_057_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00074-02', 'fileName'] <-
        'BLT00074-2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00075-04', 'fileName'] <-
        'BLT00075-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00092-06', 'fileName'] <-
        'BLT00092-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186-06', 'fileName'] <-
        'BLT00186-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186-09', 'fileName'] <-
        'BLT00186_9'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00211-03', 'fileName'] <-
        'BLT00211-3'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00211-06', 'fileName'] <-
        'BLT00211-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00211-06', 'fileName'] <-
        'BLT00211-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00214-02', 'fileName'] <-
        'BLT00214-2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00214-05', 'fileName'] <-
        'BLT00214-5'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00230-02', 'fileName'] <-
        'BLT00230-2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00243-07', 'fileName'] <-
        'BLT00243-7'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00244-04', 'fileName'] <-
        'BLT00244-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00244-06', 'fileName'] <-
        'BLT00244-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00244-06', 'fileName'] <-
        'BLT00244-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00254-05', 'fileName'] <-
        'BLT00254-5'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00254-07', 'fileName'] <-
        'BLT00254-7'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00265-04', 'fileName'] <-
        'BLT00265-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00265-04', 'fileName'] <-
        'BLT00265-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00265-08', 'fileName'] <-
        'BLT00265-8'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00271-04', 'fileName'] <-
        'BLT00271-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00271-07', 'fileName'] <-
        'BLT00271-7'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00274-06', 'fileName'] <-
        'BLT00274-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00285-03', 'fileName'] <-
        'BLT00285-3'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00285-05', 'fileName'] <-
        'BLT00285-5'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00285-05', 'fileName'] <-
        'BLT00285-5'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00286-04', 'fileName'] <-
        'BLT000286-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00286-06', 'fileName'] <-
        'BLT00286-6'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00297_02', 'fileName'] <-
        'BLT00297_2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT0244_02', 'fileName'] <-
        'BLT00244_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'QS_024-02', 'fileName'] <-
        'QS_024-2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00074-4', 'fileName'] <-
        'BLT00074_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00075-6', 'fileName'] <-
        'BLT00075_06'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186-2', 'fileName'] <-
        'BLT00186_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186-9', 'fileName'] <-
        'BLT00186_9'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00261-2', 'fileName'] <-
        'BLT00261_2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00297-2', 'fileName'] <-
        'BLT00297_2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00075-4_R1', 'fileName'] <-
        'BLT00075-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00211-3 ', 'fileName'] <-
        'BLT00211-3'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00297-2', 'fileName'] <-
        'BLT00297_2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT000214-5', 'fileName'] <-
        'BLT00214-5'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00057-4', 'fileName'] <-
        'BLT00057_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00092_4', 'fileName'] <-
        'BLT00092_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00198_2', 'fileName'] <-
        'BLT00198_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00214-4', 'fileName'] <-
        'BLT00214_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00230-4', 'fileName'] <-
        'BLT00230_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00242-2', 'fileName'] <-
        'BLT00242_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00243-5', 'fileName'] <-
        'BLT00243_05'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00274-2', 'fileName'] <-
        'BLT00274_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00243-5', 'fileName'] <-
        'BLT00243_05'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00274-2', 'fileName'] <-
        'BLT00274_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00274-4', 'fileName'] <-
        'BLT00274-05'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00286-5', 'fileName'] <-
        'BLT000286-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BAS057_02', 'fileName'] <-
        'BAS_057_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BAS033_02', 'fileName'] <-
        'BAS_033_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BAS057_02', 'fileName'] <-
        'BAS_057_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00074-4', 'fileName'] <-
        'BLT00074_04'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00075-6', 'fileName'] <-
        'BLT00075_06'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186_2', 'fileName'] <-
        'BLT00186_02'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BAS00101', 'fileName'] <-
        'BAS101'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00186-9', 'fileName'] <-
        'BLT00186_9'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00261_02', 'fileName'] <-
        'BLT00261_2'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00274-5', 'fileName'] <-
        'BLT00274-05'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00286-4', 'fileName'] <-
        'BLT000286-4'
      medianResultsIndividualClusters[medianResultsIndividualClusters$fileName ==  'BLT00297-2', 'fileName'] <-
        'BLT00297_2'

      medianResultsIndividualClusters <- medianResultsIndividualClusters[medianResultsIndividualClusters$clusterId %in% unique(dsResults$cluster_id), ]

      head(medianResultsIndividualClusters)
      medianResultsIndividualClusters$clusterIdName <-
        paste(
          medianResultsIndividualClusters$panel,
          medianResultsIndividualClusters$markersOrCell,
          medianResultsIndividualClusters$clusterName,
          medianResultsIndividualClusters$clusterId,
          sep = "_"
        )

      medianResultsIndividualClustersT <-
        reshape(
          medianResultsIndividualClusters[, c("fileName", "clusterIdName", "GPR32_median")],
          idvar = "fileName",
          timevar = "clusterIdName",
          direction = "wide"
        )

      write.csv(medianResultsIndividualClustersT, paste0("data/medianValues/",directory, clusterName, markersOrCell, "_medianValues.csv"), row.names = FALSE)
    }
  }
}
