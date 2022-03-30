try(source("R/03_functions.R"))

loadlibraries()

bCellsDirectoryName <- "bCells"
bCellsColumnNames <- c("GPR32...AF488.A","CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                       "Zombie.NIR.A","CD24...BV605.A", "CD27...BV650.A")

test <- FALSE

directoryName <- bCellsDirectoryName
columnNames <- bCellsColumnNames

preprocessing(directoryName, columnNames, test)
