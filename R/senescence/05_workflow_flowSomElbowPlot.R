try(source("R/01_functions.R"))

loadlibraries()

directoryName <- "senescence"
columnNames <- c("CD27.BV421.A", "CD45RA.BV605.A", "CD28.BV785.A", "KLRG1.PE.A",
                "CD4.PE.CF594.A", "CCR7.PE.Cy7.A", "CD8.PerCP.Cy5.5.A")
clusterName <- "meta_clusters_flowsom"
numberOfClusters <- seq(3,12)

elbowPlot(directoryName, columnNames, numberOfClusters)
