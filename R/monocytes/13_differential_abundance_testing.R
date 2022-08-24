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

markersOrCells <- c("Clusters", "CellPopulations", "Markers")

clusterNames <-
  c(
    "clusters_flowsom",
    "clusters_phenograph",
    "clusters_fast_pg",
    "meta_clusters_flowsom"
  )

n.cores <- 19
my.cluster <- parallel::makeCluster(
  n.cores
  )
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()


foreach(clusterName = clusterNames) %dopar% {
  try(source("R/01_functions.R"))

  loadlibraries()

  foreach(markersOrCell = markersOrCells) %dopar% {
    try(source("R/01_functions.R"))

    loadlibraries()
    performAllDifferentialAbundanceTests(directoryName, columnNames, clusterName, markersOrCell)
  }
}
