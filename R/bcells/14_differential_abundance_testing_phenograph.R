# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "bCells"
columnNames <-
  c(
    "fileName",
    "IgD...PerCP.Cy5.5.A",
    "CD24...BV605.A",
    "CD27...BV650.A",
    "GPR32...AF488.A",
    "FPRL1...AF647.A"
  )

clusterName <- "clusters_phenograph"

performAllDifferentialAbundanceTests(directoryName, columnNames, clusterName)
