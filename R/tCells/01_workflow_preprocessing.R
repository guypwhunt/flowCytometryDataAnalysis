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
automatedcofactors <-
  c(
    500,
    15140.92607,
    1634.74308,
    25,
    25,
    1000,
    500,
    343.22777,
    90.43229,
    2838.94867
  )

test <- FALSE

directoryName <- tCellDirectoryName
columnNames <- tCellColumnNames

preprocessing(directoryName, columnNames, test,
              automatedcofactors = automatedcofactors)
