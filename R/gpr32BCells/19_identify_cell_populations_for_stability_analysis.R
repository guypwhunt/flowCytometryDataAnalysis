try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"

markersOrCell <- "Markers"

cutoff <- gpr32BCellsCutoff

columnNames <- gpr32BCellsClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
df <- as.data.frame(df)

identifyFlowSomBoostrappedCellPopulations(df,
                                          directoryName,
                                          iterations,
                                          markersOrCell,
                                          cutoff,
                                          columnNames)
