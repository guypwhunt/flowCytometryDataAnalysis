try(source("R/01_functions.R"))

loadlibraries()

senescenceDirectoryName <- "senescence"
senescenceColumnNames <- c("CD27.BV421.A",
                           "CD45RA.BV605.A",
                           "CD28.BV785.A",
                           "KLRG1.PE.A",
                           "CD4.PE.CF594.A",
                           "CD8.BV650.A",
                           "CCR7.PE.Cy7.A",
                           "CD8.PerCP.Cy5.5.A",
                           "GPR32.AF488.A",
                           "FPRL1.AF647.A",
                           "Zombie.NIR.A")

test <- FALSE

directoryName <- senescenceDirectoryName
columnNames <- senescenceColumnNames

columnNames <- columnNames[columnNames!= "CD8.PerCP.Cy5.5.A"]
columnNames <- columnNames[columnNames!= "CD8.BV650.A"]
columnNames <- columnNames[columnNames!= "CD4.PE.CF594.A"]
columnNames <- columnNames[columnNames!= "CD45RA.BV605.A"]
columnNames <- columnNames[columnNames!= "GPR32.AF488.A"]
columnNames <- columnNames[columnNames!= "FPRL1.AF647.A"]
columnNames <- columnNames[columnNames!= "Zombie.NIR.A"]

numberOfClusters <- 24
flowsomClustering(directoryName, columnNames, numberOfClusters, test)
