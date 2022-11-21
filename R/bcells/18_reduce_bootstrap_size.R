try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "bCells"

reduceBoostrappedDataSize(directoryName, iterations)
