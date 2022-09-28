try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "monocytes"
columnNames <- monocytesColumnNames

markersOrCells <- markersOrCellsClassification

#markersOrCell <- markersOrCells[1]

clusterNames <- clusterColumns

#clusterName <- clusterNames[3]

df <-
  read.csv(paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    calculateMediansValue(directoryName, columnNames, markersOrCell,
                          clusterName, df)
  }
}
