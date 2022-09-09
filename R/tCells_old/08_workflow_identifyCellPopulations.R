try(source("R/01_functions.R"))
loadlibraries()

clusterColumns <-
  c(
    "clusters_flowsom",
    "meta_clusters_flowsom",
    "clusters_phenograph",
    "clusters_fast_pg"
  )

# tCells
columnNames <-
  c(
    "CD127.BV510.A",
    "CD8.BV650.A",
    "CD25.BV786.A",
    "FoxP3.PE.A",
    "CD45RO.PE.CF595.A",
    "CD4.PerCP.Cy5.5.A"
  )
cutoff <- c(0, 3.25, 3, 0, 2.5, 1)
directoryName <- "tCells"

#clusterName <- clusterColumns[1]
markersOrCells <- c("CellPopulations", "Markers")

for (markersOrCell in markersOrCells) {
  for (clusterColumn in clusterColumns) {
    calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
