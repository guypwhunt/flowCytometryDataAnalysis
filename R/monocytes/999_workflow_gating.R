try(source("R/01_functions.R"))

loadlibraries()

senescenceDirectoryName <- "monocytes"
test <- FALSE

directoryName <- senescenceDirectoryName
columnNames = c("CD14...BV605.A","CD16...PE.CF595.A")
cutoff = c(0,0)

gateTwoMarkersCombined(directoryName, columnNames, cutoff)
