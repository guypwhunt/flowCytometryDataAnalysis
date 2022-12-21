try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Monocytes"
columnNames <- gpr18MonocytesClusteringColumnNames

umapDimReduction(directoryName, columnNames, knn/2)
