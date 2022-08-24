# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "bCells"
columnNames <-
  c(
    "CD27...BV650.A",
    "CD24...BV605.A",
    "IgD...PerCP.Cy5.5.A"
  )

prettyColumnNames <-
  c(
    "CD27",
    "CD24",
    "IgD"
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

#clusterNames <-c(clusterNames[4])
#clusterName <- clusterNames[1]
#markersOrCellS <- c(markersOrCells[3])
#markersOrCell <- c(markersOrCells[1])


for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}

markerType <- "Functional"

columnNames <-
  c(
    "GPR32...AF488.A",
    "FPRL1...AF647.A"
  )

prettyColumnNames <-   c(
  "GPR32",
  "FPRL1"
)

for (markersOrCell in markersOrCells) {
  for (clusterName in clusterNames) {
    generateHeatmap(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames)
  }
}
