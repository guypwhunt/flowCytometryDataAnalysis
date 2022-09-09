# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "senescence"
columnNames <- c(
  "fileName",
  "CD27.BV421.A",
  "CD45RA.BV605.A",
  "CD28.BV785.A",
  "KLRG1.PE.A",
  "CD4.PE.CF594.A",
  "CD8.BV650.A",
  "CCR7.PE.Cy7.A",
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
