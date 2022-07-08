try(source("R/01_functions.R"))

loadlibraries()

bCellsDirectoryName <- "bCells"
bCellsColumnNames <- c("GPR32...AF488.A", "FPRL1...AF647.A",
                       "CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                       "Zombie.NIR.A","CD24...BV605.A", "CD27...BV650.A")

test <- FALSE

directoryName <- bCellsDirectoryName
columnNames <- bCellsColumnNames
columnNames <- columnNames[columnNames!= "Zombie.NIR.A"]
columnNames <- columnNames[columnNames!= "CD19...PE.CF595.A"]
columnNames <- columnNames[columnNames!= "GPR32...AF488.A"]
columnNames <- columnNames[columnNames!= "FPRL1...AF647.A"]

numberOfClusters <- 6
knn <- 1000
diffusionMapDimReduction(directoryName, columnNames, knn)
