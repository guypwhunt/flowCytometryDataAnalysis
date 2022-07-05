# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "monocytes"
columnNames <- c(
  "fileName",
  "CD11b...17BV421.A",
  "CD14...BV605.A",
  "CD16...PE.CF595.A",
  "CD11b.activated...PE.Cy7.A",
  "GPR32...AF488.A",
  "FPRL1...AF647.A"
)


clusterName <- "clusters_flowsom"

performAllDifferentialAbundanceTests(directoryName, columnNames, clusterName)
