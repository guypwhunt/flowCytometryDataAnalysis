try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"
columnNames <- gpr18TCellsClusteringColumnNames

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom")

numberOfClusters <- 20

generateSubsampledFlowsomClusters(directoryName, columnNames, clusterNames,
                                  numberOfClusters, iterations)
