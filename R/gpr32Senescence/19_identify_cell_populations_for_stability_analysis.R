try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32Senescence"

markersOrCell <- "Markers"

cutoff <- gpr32SenescenceCutoff

columnNames <- gpr32SenescenceClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
df <- as.data.frame(df)

identifyFlowSomBoostrappedCellPopulations(df,
                                          directoryName,
                                          iterations,
                                          markersOrCell,
                                          cutoff,
                                          columnNames)
