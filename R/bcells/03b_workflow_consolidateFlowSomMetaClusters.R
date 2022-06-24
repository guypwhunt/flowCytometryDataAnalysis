try(source("R/01_functions.R"))
loadlibraries()

directoryName <- "bCells"
columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
clusterName <- "meta_clusters_flowsom"
numberOfClusters <- seq(3,12)

consolidateFlowSomClusters(directoryName, columnNames, clusterName, numberOfClusters)
