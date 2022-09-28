try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
clusterNames <- clusterColumns

markersOrCells <- markersOrCellsClassification

# clusterName <- clusterNames[3]
# markersOrCell <- markersOrCells[1]
# panel <- "monocytes"
# fileName <- "BAS032_02"

results <- data.frame(
  clusterName = character(),
  markersOrCell = character(),
  panel = character(),
  fileName = character(),
  medianValue = double()
)

dir.create("data/medianValues", showWarnings = FALSE)

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

    for (panel in unique(dsResults$panel)) {
      df <- read.csv(paste0(
        "data/",
        panel,
        "/clusteringOutput/clusteringOutputs.csv"
      ))

      df <- df[df[, clusterName] %in% dsResults$cluster_id,]

      for (fileName in unique(df[, "fileName"])) {
        medianValue <- median(df[df[, "fileName"] == fileName, "GPR32"])

        results[nrow(results) + 1,] <-
          c(clusterName,
            markersOrCell,
            panel,
            fileName,
            medianValue)
      }
      write.csv(results, "data/medianValues/medianValues.csv", row.names = FALSE)
    }
  }
}
