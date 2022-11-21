try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "tCells"

reduceBoostrappedDataSize(directoryName, iterations)
