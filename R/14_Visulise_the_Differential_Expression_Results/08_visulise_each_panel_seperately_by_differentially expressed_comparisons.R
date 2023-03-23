try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns

markersOrCells <- markersOrCellsClassification
figureNames <-
  c(
    "DifferentialStatesStatisticscsv"#,
    # "DifferentialAbundanceStatisticscsv"
    )

markersOrCells <- markersOrCells[3]
clusterNames <- clusterNames[4]

comparisonOrder <- c(
  "Baseline pwALS vs NNC",
  "Baseline A-B vs Baseline A-L",
  "Baseline A-B vs NNC",
  "Baseline A-F vs Baseline A-S",
  "Baseline A-F vs NNC",
  "Baseline A-FB vs NNC",
  "Baseline A-FL vs NNC",
  "V2 pwALS vs Baseline pwALS",
  "Baseline A-L vs NNC",
  "Baseline A-S vs NNC",
  "Baseline A-SB vs NNC",
  "Baseline A-SL vs NNC",
  "V3 pwALS vs Baseline pwALS",
  "V3 pw vs V2 pwALS"
)

directoryNames <- c(
  # "BCells"#,
  # "Monocytes"#,
  "Senescence"#,
  # "TCells"
  )

clusterName <- "clusters_phenograph"
markersOrCell <- "Markers"
figureName <- "DifferentialStatesStatisticscsv"
directoryName <- "Senescence"

regulation <- c(
  "Upregulated"  = "Red",
  "Not Statistically Different" = "Grey",
  "Downregulated" = "Blue"
)

for (directoryName in directoryNames) {
  cellPopulationOrder <- identifyCellPopulationOrder(directoryName)

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

        yLimit <- ceiling(max(combinedDf$MinusLogFDRAdjustedPValue))

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



        combinedDf <- left_join(data.frame(Comparison = comparisonOrder), combinedDf, by = "Comparison", multiple =
                                  "all")

        combinedDf <- left_join(data.frame(typeOfCells = cellPopulationOrder), combinedDf, by = "typeOfCells", multiple =
                                  "all")

        combinedDf <- combinedDf[!is.na(combinedDf), ]

        combinedDf$typeOfCells <- gsub("\\(", "\n\\(", combinedDf$typeOfCells)

        combinedDf$regulation <-
          factor(
            combinedDf$regulation,
            levels = names(regulation)
          )

        try(combinedDf <- filterFlowCytometryComparisons(combinedDf))

        combinedDf$Panel <- factor(combinedDf$Panel)

        for (marker in na.omit(unique(combinedDf$Marker))) {
          try({
            message(marker)
            combinedDf2 <- combinedDf[combinedDf$Marker == marker,]

            combinedDf2 <-
              combinedDf2[combinedDf2$Panel %in% paste0(tolower(marker), directoryName),]

            combinedDf2 <- combinedDf2[combinedDf2$Comparison %in% unique(combinedDf2[combinedDf2$FDRAdjustedPValue< 0.05 , "Comparison"]),]

            comaprisonNames <- unique(combinedDf2$Comparison)

            numberOfPlots <- ceiling(length(comaprisonNames)/5)

            parentCells <- c("B", "M", "S", "T")

            combinedDf2$typeOfCells <-
              paste0(combinedDf2$typeOfCells, " (", as.numeric(factor(combinedDf2$typeOfCells, levels = unique(combinedDf2$typeOfCells))) - 1, ")")

            combinedDf2$typeOfCells <- factor(combinedDf2$typeOfCells, levels = unique(combinedDf2$typeOfCells))

            combinedDfGoldenSource <- combinedDf2

            for(number in seq(numberOfPlots)){
              message(number)
              combinedDf3 <- combinedDfGoldenSource

              selectedComparisons <- comaprisonNames[seq((number-1)*5 + 1, min(number * 5, length(comaprisonNames)))]

              combinedDf3 <- combinedDf3[combinedDf3$Comparison %in% selectedComparisons, ]

              comparisons <- c(
                "Baseline pwALS vs NNC" = 21,
                "Baseline A-F vs NNC" = 22,
                "Baseline A-B vs NNC" = 23,
                "Baseline A-FB vs NNC" = 24,
                "V2 pwALS vs Baseline pwALS" = 25
              )

              comparisons <- comparisons[names(comparisons) %in% selectedComparisons]

              newComparisons <- selectedComparisons[!selectedComparisons %in% names(comparisons)]

              shapesList <- c(21,22,23,24,25)

              if(length(newComparisons) > 0) {
                for(newComparison in newComparisons) {
                  value <-
                    min(shapesList[!shapesList %in% comparisons])
                  names(value) <- newComparison
                  comparisons <- c(comparisons, value)
                }
              }

              combinedDf3$Comparison <- factor(combinedDf3$Comparison, levels = names(comparisons))

            fig <-
              ggplot(
                combinedDf3,
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
                  order = 2,
                  override.aes = list(shape = 1, size = 4)
                ),
                fill = guide_colourbar(title = "log2(Fold Change)", order = 3),
                shape = guide_legend(
                  title = "Comparison",
                  order = 1,
                  override.aes = list(size = 4)
                )
              ) +
              xlab("Cell Populations") +
              ylab("-log10(Adjusted P-Value)") +
              ylim(0, yLimit) +
              scale_fill_gradientn(
                limits = c(0 - scaleLimits, scaleLimits),
                colours = c("#0000FF", "#2E2EFF", "#5C5CFF", "#8A8AFF", "#ffffff",
                            "#FF8A8A", "#FF5C5C", "#FF2E2E", "#FF0000")
              ) +
              geom_hline(yintercept = 0 - log10(0.05),
                         linetype = "dashed")

            print(fig)
            }
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
