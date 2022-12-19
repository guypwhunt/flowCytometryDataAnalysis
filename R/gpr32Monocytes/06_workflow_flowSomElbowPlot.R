try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Monocytes"
columnNames <- gpr32MonocytesClusteringColumnNames

elbowPlot(directoryName, columnNames, numberOfClusters)
