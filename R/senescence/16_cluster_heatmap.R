try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "senescence"
columnNames <- senescenceClusteringColumnNames

prettyColumnNames <- senescenceClusteringColumnNames

clusterNames <- clusterColumns

markersOrCells <- markersOrCellsClassification

markerType <- "Phenotypic"

#clusterNames <-c(clusterNames[4])
#clusterName <- clusterNames[1]
#markersOrCellS <- c(markersOrCells[3])
#markersOrCell <- c(markersOrCells[1])


for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}

markerType <- "Functional"

columnNames <-
  c(
    "GPR32"
  )

prettyColumnNames <-   c(
  "GPR32"
)

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}
