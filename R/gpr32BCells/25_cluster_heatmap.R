try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"
columnNames <- gpr32BCellsClusteringColumnNames

clusterNames <- clusterColumns[3:4]

markersOrCells <- markersOrCellsClassification[3]

markerType <- "Phenotypic"

clusterName <- clusterNames[2]
markersOrCell <- markersOrCells[1]

# df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))
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

# markerType <- "Functional"
#
# columnNames <-
#   c(
#     "GPR32"
#   )
#
# for (markersOrCell in markersOrCells) {
#   for (clusterName in clusterNames) {
#     message(markersOrCell)
#     message(clusterName)
#     try(source("R/01_functions.R"))
#     generateHeatmap(df, clusterName, directoryName, columnNames, markersOrCell, markerType)
#   }
# }

