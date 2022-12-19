try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"
columnNames <- gpr18TCellsClusteringColumnNames

clusterName <- "clusters_phenograph"

generateSubsampledPhenographClusters(directoryName,
                                     columnNames,
                                     clusterName,
                                     knn)

