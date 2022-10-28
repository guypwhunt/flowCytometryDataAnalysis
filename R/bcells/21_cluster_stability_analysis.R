try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "bCells"

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom")

markersOrCells <- markersOrCellsClassification

identifyFlowsomClusterSimilarity(directoryName,
                                 clusterNames,
                                 markersOrCells,
                                 iterations)
