# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "tCells"
columnNames <-
  c(
    "CD127.BV510.A",
    "CD8.BV650.A",
    "CD25.BV786.A",
    "FoxP3.PE.A",
    "CD45RO.PE.CF595.A",
    "CD4.PerCP.Cy5.5.A"
  )

prettyColumnNames <- c(
  "CD127",
  "CD8",
  "CD25",
  "FoxP3",
  "CD45RO",
  "CD4"
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
    "GPR32.AF488.A",
    "FPRL1.AF647.A"
  )

prettyColumnNames <- c(
  "GPR32",
  "FPRL1"
)

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}
