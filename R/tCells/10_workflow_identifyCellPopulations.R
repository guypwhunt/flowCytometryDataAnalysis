try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

#clusterName <- clusterColumns[1]

directoryName <- "tCells"
columnNames <- tCellsClusteringColumnNames

cutoff <- c(0.5,0.6,0.5,0.4,0.5,0.6)

markersOrCells <- c("CellPopulations", "Markers")

df <-
  read.csv(paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))

for (markersOrCell in markersOrCells) {
  for (clusterColumn in clusterColumns) {
    calculateClusterMarkers(df, directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
