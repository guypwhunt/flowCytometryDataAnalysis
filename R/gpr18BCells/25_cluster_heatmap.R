try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18BCells"
columnNames <- gpr18BCellsClusteringColumnNames

prettyColumnNames <- gpr18BCellsClusteringColumnNames

clusterNames <- clusterColumns[3:4]

markersOrCells <- markersOrCellsClassification[3]

markerType <- "Phenotypic"

#clusterNames <-c(clusterNames[4])
#clusterName <- clusterNames[1]
#markersOrCellS <- c(markersOrCells[3])
#markersOrCell <- c(markersOrCells[1])

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/umapDf.csv'))
df <- as.data.frame(df)

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    message(markersOrCell)
    message(clusterName)
    try(source("R/01_functions.R"))
    generateHeatmap(df, clusterName, directoryName, columnNames, markersOrCell, markerType)
  }
}
#
# markerType <- "Functional"
#
# columnNames <-
#   c(
#     "GPR18"
#   )
#
# prettyColumnNames <-   c(
#   "GPR18"
# )
#
# for (markersOrCell in markersOrCells) {
#   for (clusterName in clusterNames) {
#     message(markersOrCell)
#     message(clusterName)
#     try(source("R/01_functions.R"))
#     generateHeatmap(df, clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
#   }
# }

