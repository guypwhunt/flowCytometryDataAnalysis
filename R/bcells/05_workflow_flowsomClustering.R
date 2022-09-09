try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "bCells"
columnNames <- bCellsClusteringColumnNames

df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))

dirFCS <- paste0("./data/", directoryName, "/dataPPOutput/scaledFcs")

fcsDf <- read.twoflowdat(dir = dirFCS[1], path_CSPLR_ST = "")

for (number in numberOfClusters) {
  df <- flowsomClustering(directoryName, columnNames, number, df, fcsDf)
}

