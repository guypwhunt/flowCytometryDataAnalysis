try(source("R/03_functions.R"))

getwd()

loadlibraries()

directoryName <- "bCells"
columnNames <- c("GPR32...AF488.A","CD19...PE.CF595.A","IgD...PerCP.Cy5.5.A",
                 "CD24...BV605.A","CD27...BV650.A","Zombie.NIR.A")
test <- TRUE

preprocessing(directoryName,columnNames, test)
