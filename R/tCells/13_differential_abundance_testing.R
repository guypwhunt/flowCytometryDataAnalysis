# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "tCells"

columnNames <- c(
  "fileName",
  "CD127.BV510.A",
  "CD8.BV650.A",
  "CD25.BV786.A",
  "FoxP3.PE.A",
  "CD45RO.PE.CF595.A",
  "CD4.PerCP.Cy5.5.A",
  "GPR32.AF488.A"
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
