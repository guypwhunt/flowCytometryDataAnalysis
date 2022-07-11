try(source("R/01_functions.R"))
loadlibraries()

clusterColumns <- c("clusters_flowsom", "meta_clusters_flowsom", "clusters_phenograph", "clusters_fast_pg")

#clusterColumns <- clusterColumns[1]

# monocytes
columnNames <-
  c(
    "CD11b...17BV421.A",
    "CD11b.activated...PE.Cy7.A",
    "CD14...BV605.A",
    "CD16...PE.CF595.A"
  )
cutoff <- c(2, 0.5, 0, 0)
directoryName <- "monocytes"

markersOrCells <- c("CellPopulations", "Markers")

for (markersOrCell in markersOrCells) {
  for (clusterColumn in clusterColumns) {
    calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
