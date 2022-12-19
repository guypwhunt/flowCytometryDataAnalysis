try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"
columnNames <- c("GPR32...AF488.A", "CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                       "Zombie.NIR.A","CD24...BV605.A", "CD27...BV650.A")

prettyColumnNames <- c("GPR32", "CD19","IgD", "Zombie","CD24", "CD27")

automatedcofactors <- c(1.503293e+00, 2.339149e+04, 4.960322e-01,
                        1.371552e+03, 6.442128e+00, 1.620503e+01)

test <- FALSE

preprocessing(directoryName, columnNames, prettyColumnNames, test,
              automatedcofactors = automatedcofactors)
