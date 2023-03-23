try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Monocytes"

markersOrCell <- "Markers"

cutoff <- gpr18MonocytesCutoff

columnNames <- gpr18MonocytesClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))
df <- as.data.frame(df)
head(df)

identifyPhenographBoostrappedCellPopulations(df,
                                             directoryName,
                                             iterations,
                                             markersOrCell,
                                             cutoff,
                                             columnNames)
