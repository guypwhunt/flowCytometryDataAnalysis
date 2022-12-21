try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

#clusterName <- clusterColumns[1]

directoryName <- "gpr18BCells"
columnNames <- gpr18BCellsClusteringColumnNames

cutoff <- gpr18BCellsCutoff

markersOrCells <- c("CellPopulations", "Markers")

df <-
  fread(file=paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))
df <- as.data.frame(df)

for (markersOrCell in markersOrCells) {
  message(markersOrCell)
  for (clusterColumn in clusterColumns) {
    message(clusterColumn)
    calculateClusterMarkers(df, directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
