try(source("R/03_functions.R"))

loadlibraries()

monocytesDirectoryName <- "monocytes"
monocytesColumnNames <- c("GPR32...AF488.A", "FPRL1...AF647.A",
                          "CD11b...17BV421.A", "CD14...BV605.A",
                          "HLA.Dr...BV650.A", "CD16...PE.CF595.A",
                          "CD11b.activated...PE.Cy7.A", "Zombie.NIR.A")

test <- FALSE

directoryName <- monocytesDirectoryName
columnNames <- monocytesColumnNames

columnNames <- columnNames[columnNames!= "HLA.Dr...BV650.A"]
columnNames <- columnNames[columnNames!= "GPR32...AF488.A"]
columnNames <- columnNames[columnNames!= "FPRL1...AF647.A"]
columnNames <- columnNames[columnNames!= "Zombie.NIR.A"]

numberOfClusters <- 6
knn <- 100
fastPGClustering(directoryName, columnNames, knn)
