try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directories <- c("gpr32BCells", "gpr32Monocytes", "gpr32TCells", "gpr32Senescence")

clusterNames <- c(#"clusters_flowsom", "meta_clusters_flowsom",
  "meta_clusters_flowsomMarker",
                  #"clusters_phenograph",
  "clusters_phenographMarker")

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

    } else if (clusterName=="clusters_phenographMarker")
    {df[, clusterName] <- as.integer(rownames(df))}

    df <- df[, !colSums(is.na(df)) == nrow(df)]

    df <- df[, c(1, seq(from=ncol(df)-99, to=ncol(df)))]

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
  combinedDf[combinedDf$panel == "gpr32BCells", "panel"] <- "B Cells"
  combinedDf[combinedDf$panel == "gpr32TCells", "panel"] <- "T Cells"
  combinedDf[combinedDf$panel == "gpr32Monocytes", "panel"] <- "Monocytes"
  combinedDf[combinedDf$panel == "gpr32Senescence", "panel"] <- "Senescent T Cells"
  combinedDf$value <- combinedDf$value

  for (directory in unique(combinedDf$panel)) {
    for(cluster in unique(combinedDf[combinedDf$panel == directory, clusterName])) {
      medianValue <-
        median(combinedDf[combinedDf[, clusterName] == cluster &
                            combinedDf$panel == directory, "value"])

      if (medianValue < 0.85) {
      #  print(directory)
      #  print(clusterName)
      #  print(cluster)
      #  print(medianValue)
      }
    }
  }

  print(ggplot(combinedDf, aes_string(x=clusterName, y="value", color = "panel")) +
          geom_boxplot(show.legend = FALSE) +
          geom_hline(yintercept=c(0.85, 0.75, 0.60), linetype=rep(c("dashed","twodash","dotted"),4), color = "black") +
          ylim(0, 1) +
          scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
          facet_wrap(~panel, scales = "free") +
          theme(axis.text.x=element_blank()) +
          ylab("Jaccard Similarity Coefficient") +
          xlab(if(clusterName == "clusters_flowsom"){
            "FlowSOM Self Organising Maps"
          } else if (clusterName == "meta_clusters_flowsom"){
            "FlowSOM Meta-Clusters"
          } else if (clusterName == "clusters_phenograph"){
            "Phenograph Clusters"
          } else if (clusterName == "clusters_phenographMarker"){
            "Phenograph Cell Populations"
          } else {
            "FlowSOM Cell Populations"
          })
        )

  rm(combinedDf)
}
