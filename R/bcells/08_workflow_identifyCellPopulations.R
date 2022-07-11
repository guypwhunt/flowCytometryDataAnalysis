try(source("R/01_functions.R"))
loadlibraries()

clusterColumns <- c("clusters_flowsom", "meta_clusters_flowsom", "clusters_phenograph", "clusters_fast_pg")

#clusterName <- clusterColumns[1]

# bCells
columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
cutoff <- c(0, 0, 0)
directoryName <- "bCells"

markersOrCells <- c("CellPopulations", "Markers")

for (markersOrCell in markersOrCells) {
  for (clusterColumn in clusterColumns) {
    calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
