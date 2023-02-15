try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Monocytes"
columnNames <- gpr18MonocytesClusteringColumnNames

prettyColumnNames <- gpr18MonocytesClusteringColumnNames

clusterNames <- clusterColumns[4]

markersOrCells <- markersOrCellsClassification[3]

markerType <- "Phenotypic"

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/umapDf.csv'))
df <- as.data.frame(df)
cellPopulationOrder <- monocytesOrder

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    message(markersOrCell)
    message(clusterName)
    try(source("R/01_functions.R"))
    generateHeatmap(df, clusterName, directoryName, columnNames, markersOrCell, markerType, cellPopulationOrder)
  }
}
