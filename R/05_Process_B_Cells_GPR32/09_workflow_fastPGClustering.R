# UPDATE THE NUMBER OF CPUs in the shellScript
try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"
columnNames <- gpr32BCellsClusteringColumnNames

fastPGClustering(directoryName, columnNames, knn)
