library(data.table)

flowsomClusterBoostrappingTest <- function(directoryName) {
  try({
    df <-
      fread(
        file = paste0(
          "./data/",
          directoryName,
          '/clusteringOutput/',
          'flowsomClusterStability.csv'
        )
      )
    df <- as.data.frame(df)

    #df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
    # df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))

    results <- data.frame(
      iteration = integer(),
      numberOfClusters = integer(),
      percentageOfData = integer(),
      numberOfMetaClustersToClusters = integer()
    )

    for (x in seq(100)) {
      metaClusterName <- paste0("meta_clusters_flowsom", x)
      clusterName <- paste0("clusters_flowsom", x)

      try({
        iteration <- x
        numberOfClusters <-
          max(unique(na.omit(df[, metaClusterName])))
        percentageOfData <-
          length(na.omit(df[, metaClusterName])) / nrow(df) * 100

        clusterNumber <- sample(seq(numberOfClusters), 1)

        numberOfMetaClustersToClusters <-
          length(unique(na.omit(df[df[, clusterName] == clusterNumber, metaClusterName])))


        results[x,] <-
          c(
            iteration,
            numberOfClusters,
            percentageOfData,
            numberOfMetaClustersToClusters
          )
      })
    }

    print("")
    print(directoryName)
    print(results)
    #return(results)
  })
}

phenographClusterBoostrappingTest <- function(directoryName) {
  try({
    df <-
      fread(
        file = paste0(
          "./data/",
          directoryName,
          '/clusteringOutput/',
          'phenographClusterStability.csv'
        )
      )
    df <- as.data.frame(df)

    #df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
    # df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/phenographClusterStability.csv'))

    results <- data.frame(
      iteration = integer(),
      numberOfClusters = integer(),
      percentageOfData = integer()
    )

    for (x in seq(25)) {
      metaClusterName <- paste0("clusters_phenograph", x)

      try({
        iteration <- x
        numberOfClusters <-
          max(unique(na.omit(df[, metaClusterName])))
        percentageOfData <-
          length(na.omit(df[, metaClusterName])) / nrow(df) * 100

        clusterNumber <- sample(seq(numberOfClusters), 1)

        results[x,] <-
          c(iteration,
            numberOfClusters,
            percentageOfData)
      })
    }

    print("")
    print(directoryName)
    print(results)
  })
}

directoryNames <- c(#"bCells", "monocytes",
  "tCells", "senescence")
clusterTypes <- c("flowsom", "phenograph")

lapply(directoryNames, flowsomClusterBoostrappingTest)
lapply(directoryNames, phenographClusterBoostrappingTest)
