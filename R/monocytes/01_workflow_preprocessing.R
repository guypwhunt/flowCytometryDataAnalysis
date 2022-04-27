try(source("R/03_functions.R"))

loadlibraries()

monocytesDirectoryName <- "monocytes"
monocytesColumnNames <- c(
  "GPR32...AF488.A",
  "FPRL1...AF647.A",
  "CD11b...17BV421.A",
  "CD14...BV605.A",
  "HLA.Dr...BV650.A",
  "CD16...PE.CF595.A",
  "CD11b.activated...PE.Cy7.A",
  "Zombie.NIR.A"
)

test <- FALSE

directoryName <- monocytesDirectoryName
columnNames <- monocytesColumnNames

preprocessing(directoryName, columnNames, test)
