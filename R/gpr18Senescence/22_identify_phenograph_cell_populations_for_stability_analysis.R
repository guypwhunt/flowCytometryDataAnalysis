try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18Senescence"

markersOrCell <- "Markers"

cutoff <- gpr18SenescenceCutoff

columnNames <- gpr18SenescenceClusteringColumnNames

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))
df <- as.data.frame(df)
head(df)

identifyPhenographBoostrappedCellPopulations(df,
                                             directoryName,
                                             iterations,
                                             markersOrCell,
                                             cutoff,
                                             columnNames)
