try(source("R/03_functions.R"))

getwd()

loadlibraries()

senescenceDirectoryName <- "senescence"
senescenceColumnNames <- c("GPR32.AF488.A","KLRG1.PE.A","CD4.PE.CF594.A",
                 "CD8.PerCP.Cy5.5.A","CCR7.PE.Cy7.A", "CD28.BV785.A",
                 "Zombie.NIR.A", "CD27.BV421.A", "CD45RA.BV605.A")

bCellsDirectoryName <- "bCells"
bCellsColumnNames <- c("GPR32...AF488.A","CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                           "Zombie.NIR.A","CD24...BV605.A", "CD27...BV650.A")

tCellsDirectoryName <- "tCells"
tCellsColumnNames <- c("GPR32","CD127","CD3", "Zombie.NIR.A","CD8", "CD25", "FoxP3", "CD45RO", "CD4")

monocytesDirectoryName <- "monocytes"
monocytesColumnNames <- c("CD11b","CD14", "Zombie.NIR.A","CD16", "CD11b activated")

test <- FALSE

preprocessing(directoryName, columnNames, test)

convertToDataFrame(directoryName, columnNames)

multipleRegressionTesting(directoryName, columnNames)
