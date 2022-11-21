try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directories <- c("bCells", "monocytes", "tCells", "senescence")

directory <- "bCells"

clusterNames <- c("clusters_flowsom", "meta_clusters_flowsom", "meta_clusters_flowsomMarker")

clusterName <- clusterNames[3]

for (clusterName in clusterNames) {
  fileName <- paste0(clusterName, "Stability.csv")
  for (directory in directories) {
    #try(dev.off())
    df <-
      read.csv(
        paste0(
          "data/",
          directory,
          "/clusteringOutput/",
          fileName
        )
      )

    if(clusterName=="meta_clusters_flowsomMarker"){
      df[, clusterName] <- as.integer(rownames(df))

    }
    columnNames <- colnames(df)

    columnNames <- columnNames[columnNames != clusterName]

    longDf <- pivot_longer(df,
                           cols = columnNames,
                           names_to = 'iteration',
                           values_to = 'value')

    longDf <- as.data.frame(longDf)

    longDf[, clusterName] <- as.factor(longDf[, clusterName])

    longDf$panel <- directory

    if (exists("combinedDf")) {
      combinedDf <- rbind(longDf, combinedDf)
    } else {
      combinedDf <- longDf
    }
  }
  combinedDf[combinedDf$panel == "bCells", "panel"] <- "B Cells"
  combinedDf[combinedDf$panel == "tCells", "panel"] <- "T Cells"
  combinedDf[combinedDf$panel == "monocytes", "panel"] <- "Monocytes"
  combinedDf[combinedDf$panel == "senescence", "panel"] <- "Senescent T Cells"
  combinedDf$value <- combinedDf$value

  print(ggplot(combinedDf, aes_string(x=clusterName, y="value", color = "panel")) +
          geom_boxplot() +
          geom_hline(yintercept=0.9, linetype="dashed", color = "red") +
          ylim(0, 1) +
          scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
          facet_wrap(~panel) +
          theme(axis.text.x=element_blank()) +
          ylab("Jaccard Index / Similarity Coefficient") +
          xlab(if(clusterName == "clusters_flowsom"){
            "FlowSOM Self Organising Maps"
          } else if (clusterName == "meta_clusters_flowsom"){
            "FlowSOM Meta-Clusters"
          } else {
            "FlowSOM Marker Populations"
          })
        )

  rm(combinedDf)
}
