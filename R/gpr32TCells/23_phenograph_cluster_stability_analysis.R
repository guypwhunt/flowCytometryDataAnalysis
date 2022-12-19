try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32TCells"

clusterNames <-
  c("clusters_phenograph", "clusters_phenographMarker")
#clusterNames <- clusterNames[2]
markersOrCells <- markersOrCellsClassification

columnNames <- gpr32TCellsClusteringColumnNames

df <-
  fread(
    file = paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/phenographClusterStability.csv'
    )
  )
df <- as.data.frame(df)

results <-
  read.csv(
    paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/',
      'clusters_phenograph_Stability.csv'
    )
  )

#iterations <- 40

for (clusterName in clusterNames) {
  identifyPhenographClusterSimilarity(df,
                                      results,
                                      directoryName,
                                      clusterName,
                                      columnNames)
}
