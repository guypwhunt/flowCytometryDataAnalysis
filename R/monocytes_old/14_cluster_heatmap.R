# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "monocytes"
columnNames <-
  c(
    "CD14...BV605.A",
    "CD16...PE.CF595.A",
    "CD11b...17BV421.A",
    "CD11b.activated...PE.Cy7.A"
  )

prettyColumnNames <- c(
  "CD14",
  "CD16",
  "CD11b",
  "Activated CD11b"
)

clusterNames <-
  c(
    "clusters_flowsom",
    "clusters_phenograph",
    "clusters_fast_pg",
    "meta_clusters_flowsom"
  )

markersOrCells <- c("CellPopulations", "Markers", "Cluster")

markerType <- "Phenotypic"

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}

markerType <- "Functional"

columnNames <-
  c(
    "GPR32...AF488.A",
    "FPRL1...AF647.A",
    "HLA.Dr...BV650.A"
  )

prettyColumnNames <- c(
  "GPR32",
  "FPRL1",
  "HLA-DR"
)

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}
