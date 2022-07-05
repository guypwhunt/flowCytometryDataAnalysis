try(source("R/01_functions.R"))

loadlibraries()

directoryName <- "tCells"
columnNames <- c("CD127.BV510.A", "CD8.BV650.A", "CD25.BV786.A", "FoxP3.PE.A",
                 "CD45RO.PE.CF595.A", "CD4.PerCP.Cy5.5.A")
clusterName <- "meta_clusters_flowsom"
numberOfClusters <- seq(3,12)

elbowPlot(directoryName, columnNames, numberOfClusters)
