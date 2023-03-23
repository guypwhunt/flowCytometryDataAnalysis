try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Monocytes"
columnNames <- c(
  "CD11b - 17BV421-A",
  "CD11b - BV421-A",
  "CD11b-17BV421-A",
  "CD11b activated - PE-Cy7-A",
  "CD11b activated-PE-Cy7-A",
  "CD14 - BV605-A",
  "CD14-BV605-A",
  "CD16 - PE-CF595-A",
  "CD16-PE-CF595-A",
  "Chem23 - APC-A",
  "FPRL1 - AF647-A",
  "receptor - AF647-A",
  "GPR18 - AF488-A",
  "GPR32 - AF488-A",
  "receptor - AF488-A",
  "HLA-Dr - BV650-A",
  "HLA-Dr-BV650-A",
  "Zombie NIR-A"
)

prettyColumnNames <- c('CD11b',
                       'CD11b',
                       'CD11b',
                       'CD11b_activated',
                       'CD11b_activated',
                       'CD14',
                       'CD14',
                       'CD16',
                       'CD16',
                       'Chem23',
                       'Chem23',
                       'Chem23',
                       'GPR18',
                       'GPR18',
                       'GPR18',
                       'HLA_DR',
                       'HLA_DR',
                       'Zombie')

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
