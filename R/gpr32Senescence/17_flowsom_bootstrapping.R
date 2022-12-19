try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Senescence"

columnNames <- gpr32SenescenceClusteringColumnNames

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom")

numberOfClusters <- 15

generateSubsampledFlowsomClusters(directoryName, columnNames, clusterNames,
                                  numberOfClusters, iterations)
