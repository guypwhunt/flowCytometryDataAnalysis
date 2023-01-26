try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"
columnNames <- gpr18TCellsColumnNames

markersOrCells <- markersOrCellsClassification[3]

clusterNames <- clusterColumns[3:4]

df <-
  fread(file = paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))
df <- as.data.frame(df)

markerName <- "gpr18"

for (clusterName in clusterNames) {
  message(clusterName)
  for (markersOrCell in markersOrCells) {
    message(markersOrCell)
    calculateMediansValue(directoryName, columnNames, markersOrCell,
                          clusterName, df, markerName)
  }
  message("")
}
