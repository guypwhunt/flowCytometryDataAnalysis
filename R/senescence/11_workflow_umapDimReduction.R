try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "senescence"
columnNames <- senescenceClusteringColumnNames

umapDimReduction(directoryName, columnNames, knn/2)
