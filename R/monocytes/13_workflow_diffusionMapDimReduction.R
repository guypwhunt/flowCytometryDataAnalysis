try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "monocytes"
columnNames <- monocytesClusteringColumnNames

diffusionMapDimReduction(directoryName, columnNames, knn*20)
