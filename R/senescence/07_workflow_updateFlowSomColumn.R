try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "senescence"

df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))

df[, "meta_clusters_flowsom"] <- df[, "meta_clusters_flowsom15"]

write.csv(df, paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'), row.names = FALSE)
