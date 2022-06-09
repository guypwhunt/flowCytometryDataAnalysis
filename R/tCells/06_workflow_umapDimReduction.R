try(source("R/01_functions.R"))

loadlibraries()

tCellDirectoryName <- "tCells"
tCellColumnNames <- c("CD127.BV510.A", "CD3.BV605.A", "CD8.BV650.A",
                      "CD25.BV786.A", "FoxP3.PE.A", "CD45RO.PE.CF595.A",
                      "CD4.PerCP.Cy5.5.A", "GPR32.AF488.A", "Zombie.NIR.A",
                      "FPRL1.AF647.A")

test <- FALSE

directoryName <- tCellDirectoryName
columnNames <- tCellColumnNames

columnNames <- columnNames[columnNames!= "Zombie.NIR.A"]
columnNames <- columnNames[columnNames!= "CD19...PE.CF595.A"]
columnNames <- columnNames[columnNames!= "GPR32...AF488.A"]
columnNames <- columnNames[columnNames!= "FPRL1.AF647.A"]

numberOfClusters <- 6
knn <- 25
umapDimReduction(directoryName, columnNames, knn)
