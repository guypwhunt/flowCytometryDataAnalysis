try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "monocytes"
columnNames <- append("fileName", monocytesColumnNames)

markersOrCells <- markersOrCellsClassification

clusterNames <- clusterColumns

# clusterName <- clusterNames[1]
# markersOrCell <- markersOrCells[2]

n.cores <- 10
my.cluster <- parallel::makeCluster(
  n.cores
  )
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))

foreach(clusterName = clusterNames) %dopar% {
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))

  loadlibraries()

  foreach(markersOrCell = markersOrCells) %dopar% {
    try(source("R/01_functions.R"))
    try(source("R/00_datasets.R"))

    loadlibraries()
    performAllDifferentialAbundanceTests(df, directoryName, columnNames, clusterName, markersOrCell)
  }
}
