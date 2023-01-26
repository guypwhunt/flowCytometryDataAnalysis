try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directories <-
  c(# "BCells"#,
    # "Monocytes"#,
    "TCells"#,
    # "Senescence"
    )

markerNames <- c("gpr18", "gpr32")

clusterNames <- c(#"clusters_flowsom", "meta_clusters_flowsom",
  "meta_clusters_flowsomMarker",
  #"clusters_phenograph",
  "clusters_phenographMarker")

Stability <- c(
  "Very High" = "#00BF7D",
  "High" = "#00B0F6",
  "Moderate" = "#D89000",
  "Low" = "#F8766D"
)

clusterName <- "clusters_phenographMarker"
directory <- "gpr32BCells"
for (markerName in markerNames) {
  message("")
  message(markerName)
  for (clusterName in clusterNames) {
    message(clusterName)
    fileName <- paste0(clusterName, "Stability.csv")
    for (directory in directories) {
      message(directory)
      directory <- paste0(markerName, directory)
      #try(dev.off())
      df <-
        read.csv(paste0("data/",
                        directory,
                        "/clusteringOutput/",
                        fileName))

      df <- df[,!colSums(is.na(df)) == nrow(df)]

      df <- df[, c(1, seq(from = ncol(df) - 99, to = ncol(df)))]

      columnNames <- colnames(df)

      columnNames <- columnNames[columnNames != clusterName]

      longDf <- pivot_longer(df,
                             cols = columnNames,
                             names_to = 'iteration',
                             values_to = 'value')

      longDf <- as.data.frame(longDf)

      longDf[, clusterName] <- as.factor(longDf[, clusterName])

      longDf$panel <- directory

      combinedDf <- longDf

      combinedDf[combinedDf$panel == "gpr32BCells", "panel"] <-
        "B Cells"
      combinedDf[combinedDf$panel == "gpr32TCells", "panel"] <-
        "T Cells"
      combinedDf[combinedDf$panel == "gpr32Monocytes", "panel"] <-
        "Monocytes"
      combinedDf[combinedDf$panel == "gpr32Senescence", "panel"] <-
        "Senescent T Cells"
      combinedDf$value <- combinedDf$value

      directory <- unique(combinedDf$panel)

      combinedDf$Stability <- "Low"

      for (cluster in unique(combinedDf[, clusterName])) {
        medianValue <-
          median(combinedDf[combinedDf[, clusterName] == cluster, "value"])

        if (medianValue > 0.85) {
          combinedDf[combinedDf[, clusterName] == cluster, "Stability"] <-
            "Very High"
        } else if (medianValue > 0.75) {
          combinedDf[combinedDf[, clusterName] == cluster, "Stability"] <-
            "High"
        } else if (medianValue > 0.60) {
          combinedDf[combinedDf[, clusterName] == cluster, "Stability"] <-
            "Moderate"
        }
      }

      combinedDf$Stability <-
        factor(combinedDf$Stability, levels = names(Stability))

      combinedDf$typeOfCells <- combinedDf$panel
      combinedDf$cluster_id <- combinedDf[, clusterName]
      combinedDf <- updateMarkerNames(combinedDf)
      combinedDf$typeOfCells <-
        gsub("\\(", "\n\\(", combinedDf$typeOfCells)
      combinedDf <- combinedDf[order(combinedDf$typeOfCells),]

      combinedDf$typeOfCells <-
        paste0(combinedDf$typeOfCells, " (", as.numeric(factor(combinedDf$typeOfCells, levels = unique(combinedDf$typeOfCells))), ")")

      combinedDf$typeOfCells <- factor(combinedDf$typeOfCells, levels = unique(combinedDf$typeOfCells))


      print(
        ggplot(
          combinedDf,
          aes_string(x = "typeOfCells", y = "value", fill = "Stability")
        ) +
          geom_boxplot(show.legend = TRUE) +
          geom_hline(
            yintercept = c(0.85, 0.75, 0.60),
            linetype = rep(c("dashed", "twodash", "dotted"), length(unique(
              combinedDf$panel
            ))),
            color = "black"
          ) +
          ylim(0, 1) +
          scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
          scale_fill_manual(values = Stability, drop = FALSE) +
          theme_bw() +
          theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            size = 8
          )) +
          ylab("Jaccard Similarity Coefficient") +
          xlab(if (clusterName == "clusters_flowsom") {
            "FlowSOM Self Organising Maps"
          } else if (clusterName == "meta_clusters_flowsom") {
            "FlowSOM Meta-Clusters"
          } else if (clusterName == "clusters_phenograph") {
            "Phenograph Clusters"
          } else if (clusterName == "clusters_phenographMarker") {
            "Phenograph Cell Populations"
          } else {
            "FlowSOM Cell Populations"
          })
      )
    }
  }
}
