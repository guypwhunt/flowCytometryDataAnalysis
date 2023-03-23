try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18BCells"
columnNames <- c("GPR32 - AF488-A", "GPR18 - AF488-A",
                 #"CD19 - PE-CF595-A",
                 "IgD - PerCP-Cy5.5-A",
                 #"Zombie NIR-A",
                 "CD24 - BV605-A", "CD27 - BV650-A")

prettyColumnNames <- c("GPR18", "GPR18",
                       #"CD19",
                       "IgD",
                       #"Zombie",
                       "CD24", "CD27")

automatedcofactors <- c(2.166595e+01,
                        #1.717232e+04,
                        3.230456e-01,
                        #1.407461e+03,
                        1.928897e+00,
                        2.061957e+00)

test <- FALSE

preprocessing(directoryName, columnNames, prettyColumnNames, test,
              automatedcofactors = automatedcofactors)
