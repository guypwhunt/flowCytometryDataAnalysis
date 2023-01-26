try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom",
  "meta_clusters_flowsomMarker")

markersOrCells <- markersOrCellsClassification

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
df <- as.data.frame(df)

results <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/', 'meta_clusters_flowsom_Stability.csv'))

# clusterName <- clusterNames[3]

for (clusterName in clusterNames) {
  identifyFlowsomClusterSimilarity(df,
                                   results,
                                   directoryName,
                                   clusterName,
                                   iterations)
}
