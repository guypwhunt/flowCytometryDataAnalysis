try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Monocytes"
columnNames <- c(
  "GPR32...AF488.A",
  "CD11b...17BV421.A",
  "CD14...BV605.A",
  "HLA.Dr...BV650.A",
  "CD16...PE.CF595.A",
  "CD11b.activated...PE.Cy7.A",
  "Zombie.NIR.A"
)

prettyColumnNames <- c("GPR32", "CD11b","CD14", "HLA_DR","CD16", "CD11b_activated", "Zombie")

automatedcofactors <- c(130.957690, 28828.722419, 41.415505,
                        1194.595206, 10.386289, 32.117761, 932.989806)

test <- FALSE
gate <- TRUE
gateTogether <- TRUE

gateColumns <- data.frame (CD14  = c(0),
                           CD16 = c(0))

preprocessing(directoryName, columnNames, prettyColumnNames,
              test, gate, gateTogether, gateColumns,
              automatedcofactors = NULL)
