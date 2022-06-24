try(source("R/01_functions.R"))
loadlibraries()

clusterColumns <- c("clusters_flowsom", "meta_clusters_flowsom")

for (clusterColumn in clusterColumns) {
  # bCells
  columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
  cutoff <- c(0, 0, 0)
  directoryName <- "bCells"
  calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff)

  # monocytes
  columnNames <- c("CD11b...17BV421.A", "CD11b.activated...PE.Cy7.A", "CD14...BV605.A", "CD16...PE.CF595.A" )
  cutoff <- c(2, 0.5, 0, 0)
  directoryName <- "monocytes"
  calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff)


  # tCells
  columnNames <- c("CD127.BV510.A", "CD8.BV650.A", "CD25.BV786.A", "FoxP3.PE.A", "CD45RO.PE.CF595.A", "CD4.PerCP.Cy5.5.A")
  cutoff <- c(0, 3.25, 3, 0, 2.5, 1)
  directoryName <- "tCells"
  calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff)

  # senescentTCells
  columnNames <- c("CD27.BV421.A", "CD45RA.BV605.A", "CD28.BV785.A", "KLRG1.PE.A", "CD4.PE.CF594.A", "CCR7.PE.Cy7.A", "CD8.PerCP.Cy5.5.A")
  cutoff <- c(2.5, 2, 0, 0, 1.5, 0, 1.75)
  directoryName <- "senescence"
  calculateClusterMarkers(directoryName, clusterColumn, columnNames, cutoff)
}
