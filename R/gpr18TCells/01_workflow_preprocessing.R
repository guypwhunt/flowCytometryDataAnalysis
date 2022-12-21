try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"

columnNames <- c(
  'CD127 - BV510-A',
  'CD127-BV510-A',
  'CD25 - BV786-A',
  'CD25-BV786-A',
  #'CD3 - BV605-A',
  #'CD3-BV605-A',
  'CD4 - PerCP-Cy5.5-A',
  'CD4-PerCP-Cy5.5-A',
  'CD45RO - PE-CF595-A',
  'CD45RO-PE-CF595-A',
  'CD8 - BV650-A',
  'CD8-BV650-A',
  'FoxP3 - PE-A',
  'FoxP3-PE-A',
  'GPR18 - AF488-A',
  'GPR32-AF488-A',
  'Zombie NIR-A'
)

prettyColumnNames <-
  c(
    'CD127',
    'CD127',
    'CD25',
    'CD25',
    #'CD3',
    #'CD3',
    'CD4',
    'CD4',
    'CD45RO',
    'CD45RO',
    'CD8',
    'CD8',
    'FoxP3',
    'FoxP3',
    'GPR18',
    'GPR18',
    'Zombie'
  )

automatedcofactors <- c(2962.46160, 3709.97756,
                        123.13927, 2536.12924, 3602.32826,
                        56.13862, 402.84477, 101.70538)

test <- FALSE

preprocessing(directoryName, columnNames, prettyColumnNames,
              test, automatedcofactors = automatedcofactors)
