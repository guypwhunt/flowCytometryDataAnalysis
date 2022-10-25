try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

#clusterName <- clusterColumns[1]

directoryName <- "bCells"

markersOrCells <- markersOrCellsClassification

#markersOrCell <- markersOrCells[1]

clusterNames <- clusterColumns

#clusterName <- clusterNames[3]

df <-
  read.csv(paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))

my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

foreach(clusterName = clusterNames, markersOrCell = markersOrCells) %dopar% {
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))

  calculateCounts(directoryName, markersOrCell,
                        clusterName, df)
}
