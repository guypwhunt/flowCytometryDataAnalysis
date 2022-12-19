try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

#clusterName <- clusterColumns[1]

directoryName <- "gpr32Senescence"
columnNames <- gpr32SenescenceColumnNames

markersOrCells <- markersOrCellsClassification

#markersOrCell <- markersOrCells[1]

clusterNames <- clusterColumns

#clusterName <- clusterNames[3]

df <-
  fread(file=paste0(
    "data/",
    directoryName,
    "/clusteringOutput/clusteringOutputs.csv"
  ))
df <- as.data.frame(df)

my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

foreach(clusterName = clusterNames, markersOrCell = markersOrCells) %dopar% {
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))

  calculateMediansValue(directoryName, columnNames, markersOrCell,
                        clusterName, df)
}
