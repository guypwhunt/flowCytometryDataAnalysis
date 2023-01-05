try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"
columnNames <- append("fileName", gpr18TCellsColumnNames)

markersOrCells <- markersOrCellsClassification[3]

clusterNames <- clusterColumns[3:4]

n.cores <- 1

# clusterName <- "meta_clusters_flowsom"
# markersOrCell <- markersOrCells[1]
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

df <-
  fread(file = paste0(
    "./data/",
    directoryName,
    '/clusteringOutput/clusteringOutputs.csv'
  ))
df <- as.data.frame(df)

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    try(source("R/01_functions.R"))
    try(source("R/00_datasets.R"))

    loadlibraries()

    # if (markersOrCell == "Markers") {
    #   if (clusterName == "meta_clusters_flowsom") {
    #     cutoffDf <- fread(
    #       file =
    #         paste0(
    #           "data/",
    #           directoryName,
    #           "/clusteringOutput/meta_clusters_flowsomMarkerStability.csv"
    #         )
    #     )
    #
    #     cutoffDf <- as.data.frame(cutoffDf)
    #
    #     cutoffDf$rowMedian <- apply(cutoffDf[,-1], 1, median)
    #
    #     colnames(cutoffDf)[1] <- clusterName
    #
    #     cutoffDf <- cutoffDf[, c(clusterName, "rowMedian")]
    #
    #     cutoffDf <<-
    #       cutoffDf[cutoffDf$rowMedian < clusterStabilityCutoff, ]
    #   } else if (clusterName == "clusters_phenograph") {
    #     cutoffDf <- fread(
    #       file =
    #         paste0(
    #           "data/",
    #           directoryName,
    #           #"/clusteringOutput/meta_clusters_flowsomMarkerStability.csv"
    #           "/clusteringOutput/clusters_phenographMarkerStability.csv"
    #           )
    #     )
    #     cutoffDf <- as.data.frame(cutoffDf)
    #
    #     cutoffDf$rowMedian <- apply(cutoffDf[,-1], 1, median)
    #
    #     colnames(cutoffDf)[1] <- clusterName
    #
    #     cutoffDf <- cutoffDf[, c(clusterName, "rowMedian")]
    #
    #     cutoffDf <-
    #       cutoffDf[cutoffDf$rowMedian < clusterStabilityCutoff, ]
    #
    #   }
    # } else if (markersOrCell != "Markers") {
    #   cutoffDf <- fread(
    #     file =
    #       paste0(
    #         "data/",
    #         directoryName,
    #         "/clusteringOutput/",
    #         clusterName,
    #         markersOrCell,
    #         "CountsOverCutoff.csv"
    #       )
    #   )
    #
    #   cutoffDf <- as.data.frame(cutoffDf)
    # }

    performAllDifferentialAbundanceTests(df,
                                         directoryName,
                                         columnNames,
                                         clusterName,
                                         markersOrCell)
  }
}
