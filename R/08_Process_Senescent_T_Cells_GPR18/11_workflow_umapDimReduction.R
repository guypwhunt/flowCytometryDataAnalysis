try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Senescence"
columnNames <- gpr18SenescenceClusteringColumnNames

umapDimReduction(directoryName, columnNames, knn/2)
