# UPDATE THE NUMBER OF CPUs in the shellScript
try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Senescence"
columnNames <- gpr32SenescenceClusteringColumnNames

fastPGClustering(directoryName, columnNames, knn)
