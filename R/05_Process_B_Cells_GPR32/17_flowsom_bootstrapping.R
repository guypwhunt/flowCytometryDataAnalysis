try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"

columnNames <- gpr32BCellsClusteringColumnNames

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom")

numberOfClusters <- 8

generateSubsampledFlowsomClusters(directoryName, columnNames, clusterNames,
                                  numberOfClusters, iterations)
