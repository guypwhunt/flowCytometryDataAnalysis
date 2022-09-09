try(source("R/01_functions.R"))

loadlibraries()

tCellDirectoryName <- "tCells"
tCellColumnNames <- c(
  "CD127.BV510.A",
  "CD3.BV605.A",
  "CD8.BV650.A",
  "CD25.BV786.A",
  "FoxP3.PE.A",
  "CD45RO.PE.CF595.A",
  "CD4.PerCP.Cy5.5.A",
  "GPR32.AF488.A",
  "Zombie.NIR.A",
  "FPRL1.AF647.A"
)


prettyColumnNames <- c("CD127", "CD3","CD8", "CD25","FoxP3",
                       "CD45RO", "CD4", "GPR32", "Zombie")

automatedcofactors_updated <-
  c(
    500,
    15140.92607,
    1634.74308,
    25,
    25,
    1000,
    500,
    343.22777,
    90.43229
  )

automatedcofactors <- c(3161.24533, 15140.92607,  1634.74308,
                        605.70310, 49.16272, 2789.74216,
                        2071.21481, 343.22777, 90.43229)


test <- FALSE

directoryName <- tCellDirectoryName
columnNames <- tCellColumnNames

preprocessing(directoryName, columnNames,
              prettyColumnNames, test,
              automatedcofactors = NULL)
