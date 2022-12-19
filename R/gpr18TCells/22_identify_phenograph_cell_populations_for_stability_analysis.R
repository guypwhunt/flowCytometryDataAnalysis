try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"
columnNames <- gpr18TCellsClusteringColumnNames

markersOrCell <- "Markers"

cutoff <- gpr18TCellsCutoff

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))
df <- as.data.frame(df)
head(df)

identifyPhenographBoostrappedCellPopulations(df,
                                             directoryName,
                                             iterations,
                                             markersOrCell,
                                             cutoff,
                                             columnNames)
