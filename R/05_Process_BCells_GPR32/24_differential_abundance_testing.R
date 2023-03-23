try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"
columnNames <- append("fileName", gpr32BCellsColumnNames)

markersOrCells <- markersOrCellsClassification[3]
message(markersOrCells)
clusterNames <-  clusterColumns[3:4]
message(clusterNames)

n.cores <- 1

my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

df <-
  fread(file=paste0(
    "./data/",
    directoryName,
    '/clusteringOutput/clusteringOutputs.csv'
  ))
df <- as.data.frame(df)

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    try(source("R/01_functions.R"))
    try(source("R/00_datasets.R"))

    loadlibraries()

    performAllDifferentialAbundanceTests(df,
                                         directoryName,
                                         columnNames,
                                         clusterName,
                                         markersOrCell)
  }
}
