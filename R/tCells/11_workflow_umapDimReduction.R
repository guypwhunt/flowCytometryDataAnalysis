try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "tCells"
columnNames <- tCellsClusteringColumnNames

umapDimReduction(directoryName, columnNames, knn/2)
