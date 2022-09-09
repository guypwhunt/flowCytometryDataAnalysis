try(source("R/01_functions.R"))

loadlibraries()

monocytesDirectoryName <- "monocytes"
monocytesColumnNames <- c(
  "GPR32...AF488.A",
  "CD11b...17BV421.A",
  "CD14...BV605.A",
  "HLA.Dr...BV650.A",
  "CD16...PE.CF595.A",
  "CD11b.activated...PE.Cy7.A",
  "Zombie.NIR.A"
)

automatedcofactors <- c(130.957690, 28828.722419, 41.415505,
                        1194.595206, 10.386289, 32.117761, 932.989806)

test <- FALSE

directoryName <- monocytesDirectoryName
columnNames <- monocytesColumnNames
gateColumns <- data.frame (CD14...BV605.A  = c(0),
                           CD16...PE.CF595.A = c(0))
gate <- TRUE
gateTogether <- TRUE

preprocessing(directoryName, columnNames, test, gate, gateTogether,
              gateColumns, automatedcofactors = NULL)
