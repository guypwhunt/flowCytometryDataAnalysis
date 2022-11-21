try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "tCells"

markersOrCell <- "Markers"

cutoff <- tCellsCutoff

columnNames <- tCellsClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
df <- as.data.frame(df)

identifyFlowSomBoostrappedCellPopulations(df,
                                          directoryName,
                                          iterations,
                                          markersOrCell,
                                          cutoff,
                                          columnNames)
