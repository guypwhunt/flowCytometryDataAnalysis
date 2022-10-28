try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "bCells"
columnNames <- append("fileName", bCellsColumnNames)

markersOrCells <- markersOrCellsClassification

clusterNames <-  clusterColumns

# clusterName <- clusterNames[1]
# markersOrCell <- markersOrCells[1]
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

df <-
  read.csv(paste0(
    "./data/",
    directoryName,
    '/clusteringOutput/clusteringOutputs.csv'
  ))

foreach(clusterName = clusterNames) %:%
  foreach(markersOrCell = markersOrCells) %dopar% {
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))

  loadlibraries()

  cutoffDf <- read.csv(
    paste0(
      "data/",
      directoryName,
      "/clusteringOutput/",
      clusterName,
      markersOrCell,
      "CountsOverCutoff.csv"
    )
  )

  performAllDifferentialAbundanceTests(df, cutoffDf, directoryName, columnNames, clusterName, markersOrCell)
}
