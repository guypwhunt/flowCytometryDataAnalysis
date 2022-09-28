try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "senescence"
columnNames <- c("CD27.BV421.A",
                 "CD45RA.BV605.A",
                 "CD28.BV785.A",
                 "KLRG1.PE.A",
                 "CD4.PE.CF594.A",
                 "CD8.PerCP.Cy5.5.A",
                 "CCR7.PE.Cy7.A",
                 "GPR32.AF488.A",
                 "Zombie.NIR.A"
)

prettyColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4", "CD8",
                       "CCR7", "GPR32", "Zombie")

automatedcofactors <- c(5975.868290, 4937.436002, 41.415505, 0.361972, 8742.139346, 1500,
                        25.120859, 45.300406, 105.314192)

test <- FALSE
gateColumns <- data.frame (CD8 = c(2),
                           CD4 = c(1.5))
gate <- TRUE
gateTogether <- TRUE

preprocessing(directoryName, columnNames, prettyColumnNames,
              test, gate, gateTogether, gateColumns,
              automatedcofactors = automatedcofactors)
