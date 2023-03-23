try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Senescence"
columnNames <- c(
  "CD27-BV421-A",
  "CD45RA-BV605-A",
  "CD28-BV785-A",
  "KLRG1-PE-A",
  "CD8-PerCP-Cy5.5-A",
  "CCR7-PE-Cy7-A",
  "GPR32-AF488-A",
  "Zombie NIR-A",
  "CD4-PE-CF594-A"
)

prettyColumnNames <-
  c(
    "CD27",
    "CD45RA",
    "CD28",
    "KLRG1",
    "CD8",
    "CCR7",
    "GPR18",
    "Zombie",
    "CD4"
    )

automatedcofactors <-
  c(
    6317.344,
    4937.436,
    79.87204,
    0.9455981,
    7641.22212,
    126.26644,
    46.45083,
    112.57901,
    8742.139346 #85.20,
  )

test <- FALSE
gateColumns <- data.frame (CD8 = c(0.75),
                           CD4 = c(1.5))
gate <- TRUE
gateTogether <- TRUE

preprocessing(directoryName, columnNames, prettyColumnNames,
              test, gate, gateTogether, gateColumns,
              automatedcofactors = automatedcofactors)
