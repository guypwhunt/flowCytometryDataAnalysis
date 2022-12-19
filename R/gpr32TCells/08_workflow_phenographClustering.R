try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32TCells"
columnNames <- gpr32TCellsClusteringColumnNames

phenographClustering(directoryName, columnNames, knn)
