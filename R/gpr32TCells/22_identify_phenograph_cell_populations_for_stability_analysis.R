try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32TCells"

markersOrCell <- "Markers"

cutoff <- gpr32TCellsCutoff

columnNames <- gpr32TCellsClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))
df <- as.data.frame(df)
head(df)

identifyPhenographBoostrappedCellPopulations(df,
                                             directoryName,
                                             iterations,
                                             markersOrCell,
                                             cutoff,
                                             columnNames)
