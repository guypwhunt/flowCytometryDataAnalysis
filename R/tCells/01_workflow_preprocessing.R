try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "tCells"
columnNames <- c(
  "CD127.BV510.A",
  "CD3.BV605.A",
  "CD8.BV650.A",
  "CD25.BV786.A",
  "FoxP3.PE.A",
  "CD45RO.PE.CF595.A",
  "CD4.PerCP.Cy5.5.A",
  "GPR32.AF488.A",
  "Zombie.NIR.A"
)

prettyColumnNames <- c("CD127", "CD3","CD8",
                       "CD25", "FoxP3","CD45RO",
                       "CD4", "GPR32", "Zombie")

automatedcofactors <- c(3161.24533, 15140.92607,  1634.74308,
                        605.70310, 49.16272, 2789.74216,
                        2071.21481, 343.22777, 90.43229)

test <- FALSE

preprocessing(directoryName, columnNames, prettyColumnNames,
              test, automatedcofactors = automatedcofactors)
