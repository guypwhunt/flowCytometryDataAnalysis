# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "senescence"
columnNames <-
  c(
    "CD4.PE.CF594.A",
    "CD8.BV650.A",
    "CD45RA.BV605.A",
    "CD27.BV421.A",
    "CD28.BV785.A",
    "KLRG1.PE.A",
    "CCR7.PE.Cy7.A"
  )

prettyColumnNames <-   c(
  "CD4",
  "CD8",
  "CD45RA",
  "CD27",
  "CD28",
  "KLRG1",
  "CCR7"
)

clusterNames <-
  c(
    "clusters_flowsom",
    "clusters_phenograph",
    "clusters_fast_pg",
    "meta_clusters_flowsom"
  )

markersOrCells <- c("CellPopulations", "Markers")

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
