try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns

markersOrCells <- markersOrCellsClassification
figureNames <-
  c("DifferentialStatesStatisticscsv"#,
    #"DifferentialAbundanceStatisticscsv"
    )

markersOrCells <- markersOrCells[3]
clusterNames <- clusterNames[4]

directoryNames <- c(
  # "BCells"#,
  "Monocytes"#,
  # "Senescence",
  # "TCells"
  )

clusterName <- "clusters_phenograph"
markersOrCell <- "Markers"
figureName <- "DifferentialStatesStatisticscsv"
directoryName <- "TCells"

comparisons <- c(
  "ALS vs Control" = 21,
  "Fast vs Controls"  = 22,
  "Fast vs Slow"  = 23,
  "Slow vs Control"  = 24,
  "Bulbar vs Control" = 25
)

regulation <- c(
  "Upregulated"  = "Red",
  "Not Statistically Different" = "Grey",
  "Downregulated" = "Blue"
)

for (directoryName in directoryNames) {
  for (clusterName in clusterNames) {
    for (markersOrCell in markersOrCells) {
      for (figureName in figureNames) {
        message(directoryName)
        message(clusterName)
        message(markersOrCell)
        message(figureName)

        combinedDf <-
          fread(
            paste0(
              "data/pValueAdjustmentsResults/",
              clusterName,
              markersOrCell,
              figureName,
              ".csv"
            )
          )
        combinedDf <- as.data.frame(combinedDf)

        colnames(combinedDf) <-
          gsub(" ", "", colnames(combinedDf), fixed = TRUE)
        colnames(combinedDf) <-
          gsub("-", "", colnames(combinedDf), fixed = TRUE)
        colnames(combinedDf) <-
          gsub("(", "", colnames(combinedDf), fixed = TRUE)
        colnames(combinedDf) <-
          gsub(")", "", colnames(combinedDf), fixed = TRUE)

        colnames(combinedDf)[colnames(combinedDf) %in% c("ID", "CellPopulationName")] <-
          c("cluster_id", "typeOfCells")


        combinedDf <- combinedDf[!is.na(combinedDf$LogFoldChange), ]


        scaleLimits <-
          ceiling(max(abs(min(
            na.omit(combinedDf$LogFoldChange)
          )), max(
            na.omit(combinedDf$LogFoldChange)
          )))

        combinedDf$regulation <- "Not Statistically Different"

        try({combinedDf[combinedDf$FDRAdjustedPValue < 0.05 &
                         combinedDf$LogFoldChange < 0, "regulation"] <-
              "Downregulated"})
        try({combinedDf[combinedDf$FDRAdjustedPValue < 0.05 &
                         combinedDf$LogFoldChange > 0, "regulation"] <-
              "Upregulated"})

        combinedDf <- updateMarkerNames(combinedDf)

        combinedDf$typeOfCells <- gsub("\\(", "\n\\(", combinedDf$typeOfCells)

        combinedDf <- combinedDf[order(combinedDf$typeOfCells),]

        combinedDf$regulation <-
          factor(
            combinedDf$regulation,
            levels = names(regulation)
          )

        combinedDf <-
          combinedDf[combinedDf$Comparison %in% names(comparisons),]

        combinedDf[combinedDf$JaccardSimilarityCoefficient < 0.65, "jaccard"] <-
          "Low"
        combinedDf[combinedDf$JaccardSimilarityCoefficient > 0.65, "jaccard"] <-
          "Moderate"
        combinedDf[combinedDf$JaccardSimilarityCoefficient > 0.75, "jaccard"] <-
          "High"
        combinedDf[combinedDf$JaccardSimilarityCoefficient > 0.85, "jaccard"] <-
          "Very High"

        combinedDf$Panel <- factor(combinedDf$Panel)
        combinedDf <-
          combinedDf[order(combinedDf$JaccardSimilarityCoefficient), ]
        combinedDf$jaccard <-
          factor(combinedDf$jaccard,
                 levels = c("Very High", "High", "Moderate", "Low"))
        combinedDf <-
          combinedDf[order(combinedDf$typeOfCells), ]
        combinedDf <-
          combinedDf[order(gsub(
            "[0123456789]",
            "",
            gsub("gpr", "", combinedDf$Panel, fixed = TRUE)
          )), ]

        for (marker in unique(combinedDf$Marker)) {
          try({
            message(marker)
            combinedDf2 <- combinedDf[combinedDf$Marker == marker,]

            combinedDf2 <-
              combinedDf2[combinedDf2$Panel %in% paste0(tolower(marker), directoryName),]

            combinedDf2 <- combinedDf2[order(combinedDf2$typeOfCells), ]

            parentCells <- c("B", "M", "S", "T")

            typeOfCellsFinalLevel <- parentCells[parentCells %in% unique(combinedDf2$typeOfCells)]
            typeOfCellsLevels <- unique(combinedDf2$typeOfCells)[!unique(combinedDf2$typeOfCells) %in% c("B", "M", "S", "T")]

            typeOfCellsLevels <- c(typeOfCellsLevels, typeOfCellsFinalLevel)

            combinedDf2 <- left_join(data.frame(typeOfCells = typeOfCellsLevels), combinedDf2, by = "typeOfCells")

            combinedDf2$typeOfCells <-
              paste0(combinedDf2$typeOfCells, " (", as.numeric(factor(combinedDf2$typeOfCells, levels = typeOfCellsLevels)), ")")

            combinedDf2$typeOfCells <- factor(combinedDf2$typeOfCells, levels = unique(combinedDf2$typeOfCells))

            fig <-
              ggplot(
                combinedDf2,
                aes(
                  x = as.factor(typeOfCells),
                  y = MinusLogFDRAdjustedPValue,
                  fill = LogFoldChange,
                  colour = regulation,
                  shape = Comparison,
                  size = jaccard
                )
              ) +
              geom_point(
                alpha = 1,
                size = 5,
                stroke = 1,
                position = position_jitter(width = 0.1, height = 0)
              ) +
              theme_bw() +
              theme(
                axis.text.x = element_text(
                  angle = 90,
                  vjust = 0.5,
                  hjust = 1,
                  size = 8
                ),
                legend.title = element_text(size = 9),
                legend.text = element_text(size = 8),
                legend.justification = "top"
              ) +
              scale_shape_manual(values = comparisons) +
              scale_colour_manual(values = regulation, drop = FALSE) +
              guides(
                color = guide_legend(
                  title = if (figureName == "DifferentialStatesStatisticscsv") {
                    "Expression"
                  } else {
                    "Counts"
                  }
                  ,
                  order = 1,
                  override.aes = list(shape = 1, size = 4)
                ),
                fill = guide_colourbar(title = "log2(Fold Change)", order = 2),
                shape = guide_legend(
                  title = "Comparison",
                  order = 3,
                  override.aes = list(size = 4)
                )
              ) +
              xlab("Cell Populations") +
              ylab("-log10(Adjusted P-Value)") +
              ylim(0, ceiling(max(combinedDf$MinusLogFDRAdjustedPValue))) +
              geom_hline(yintercept = 0 - log10(0.05),
                         linetype = "dashed") +
              scale_fill_gradientn(
                limits = c(0 - scaleLimits, scaleLimits),
                colours = c("#0000FF", "#2E2EFF", "#5C5CFF", "#8A8AFF", "#ffffff",
                            "#FF8A8A", "#FF5C5C", "#FF2E2E", "#FF0000")
              )

            print(fig)
          })
        }

        message("")

        rm(list = ls()[!ls() %in% c(
          "clusterNames",
          "markersOrCells",
          "figureNames",
          "directoryNames",
          "clusterName",
          "markersOrCell",
          "figureName",
          "directoryName",
          "comparisons",
          "regulation"
        )])
        try(source("R/01_functions.R"))
        try(source("R/00_datasets.R"))
      }
    }
  }
}
