try(source("R/03_functions.R"))

loadlibraries()

bCellsDirectoryName <- "bCells"
bCellsColumnNames <- c("GPR32...AF488.A","CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                       "Zombie.NIR.A","CD24...BV605.A", "CD27...BV650.A")

test <- FALSE

directoryName <- bCellsDirectoryName
columnNames <- bCellsColumnNames
columnNames <- columnNames[columnNames!= "Zombie.NIR.A"]
columnNames <- columnNames[columnNames!= "CD19...PE.CF595.A"]
columnNames <- columnNames[columnNames!= "GPR32...AF488.A"]

numberOfClusters <- 6
knn <- 50
phenographClustering(directoryName, columnNames, knn)
