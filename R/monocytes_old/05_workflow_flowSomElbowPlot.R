try(source("R/01_functions.R"))

loadlibraries()

directoryName <- "monocytes"
columnNames <- c("CD11b...17BV421.A", "CD14...BV605.A",
                 "CD16...PE.CF595.A", "CD11b.activated...PE.Cy7.A")
clusterName <- "meta_clusters_flowsom"
numberOfClusters <- seq(3,12)

elbowPlot(directoryName, columnNames, numberOfClusters)
