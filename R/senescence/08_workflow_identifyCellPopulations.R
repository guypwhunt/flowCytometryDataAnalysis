try(source("R/01_functions.R"))
loadlibraries()

clusterColumns <- c("clusters_flowsom", "meta_clusters_flowsom", "clusters_phenograph", "clusters_fast_pg")

# senescentTCells
columnNames <-
  c(
    "CD27.BV421.A",
    "CD45RA.BV605.A",
    "CD28.BV785.A",
    "KLRG1.PE.A",
    "CD4.PE.CF594.A",
    "CCR7.PE.Cy7.A",
    "CD8.PerCP.Cy5.5.A"
  )
cutoff <- c(2.5, 2, 0, 0, 1.5, 0, 1.75)
directoryName <- "senescence"

markersOrCells <- c("CellPopulations", "Markers")

for (markersOrCell in markersOrCells) {
  for (clusterColumn in clusterColumns) {
    calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff, markersOrCell)
  }
}
