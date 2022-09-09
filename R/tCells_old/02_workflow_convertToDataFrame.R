try(source("R/01_functions.R"))

loadlibraries()

tCellDirectoryName <- "tCells"
tCellColumnNames <- c("CD127","CD8", "CD25","FoxP3",
                      "CD45RO", "CD4", "GPR32")

columnsToScale <- tCellColumnNames[tCellColumnNames != "GPR32"]

test <- FALSE

directoryName <- tCellDirectoryName
columnNames <- tCellColumnNames

convertToDataFrame(directoryName, columnNames, columnsToScale, test)
