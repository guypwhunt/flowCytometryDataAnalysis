try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Monocytes"

columnNames <- gpr32MonocytesClusteringColumnNames

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom")

numberOfClusters <- 15

generateSubsampledFlowsomClusters(directoryName, columnNames, clusterNames,
                                  numberOfClusters, iterations)
