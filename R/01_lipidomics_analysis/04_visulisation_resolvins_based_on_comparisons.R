try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

figureDirectory <- "data/lipidomics/Results/"

fileNames <- list.files(figureDirectory, pattern =".csv")

fileNamesPath <- paste0(figureDirectory, fileNames)

dfs <- lapply(fileNamesPath, read.csv)

names(dfs) <- fileNames

i <- 1

for(df in dfs){
  df$experiment <- names(dfs)[i]
  i <- i + 1

  if (exists("combinedDf")){
    combinedDf <- rbind(combinedDf, df)
  } else {
    combinedDf <- df
  }
}

summary(combinedDf$logFC)

df <- combinedDf

combinedDf$lipid <- factor(combinedDf$X, levels = unique(combinedDf[order(combinedDf$X), "X"]))

combinedDf <- filterLipidomicsComparisons(combinedDf)

combinedDf$regulation <- "Not Statistically Different"

try({combinedDf[combinedDf$adj.P.Val < 0.05 &
                  combinedDf$logFC < 0, "regulation"] <-
  "Downregulated"})
try({combinedDf[combinedDf$adj.P.Val < 0.05 &
                  combinedDf$logFC > 0, "regulation"] <-
  "Upregulated"})

scaleLimits <-
  ceiling(max(abs(min(
    na.omit(df$logFC)
  )), max(
    na.omit(df$logFC)
  )))

regulation <- c(
  "Upregulated"  = "Red",
  "Not Statistically Different" = "Grey",
  "Downregulated" = "Blue"
)

combinedDf$regulation <- factor(combinedDf$regulation, levels = names(regulation))

comparisonsList <- list(
  c(
    "V2 A-F vs V2 A-S" = 21,
    "V2 A-FB vs V2 A-SB"  = 22,
    "V2 A-FL vs V2 A-SL" = 23,
    "V2 A-SB vs Baseline A-SB"  = 24
  )

  # alsAndSubGroupsVsControlsComparisons,
  # progressionSubGroupsVsControlsComparisons,
  # progressionSubGroupsVsSubGroupsComparisons,
  # londitudinalComparisons[1]
  )

combinedDfGoldenSource <- combinedDf

for(comparisons in comparisonsList) {
  combinedDf <- combinedDfGoldenSource
  combinedDf <-
    combinedDf[combinedDf$experiment %in% names(comparisons),]
  combinedDf <- combinedDf[grep("Rv", combinedDf$X), ]
  combinedDf <- combinedDf[nchar(combinedDf$X) == 4, ]
  combinedDf$experiment <-
    factor(combinedDf$experiment, levels = names(comparisons))

  print(
    ggplot(
      combinedDf,
      aes(
        x = lipid,
        y = -log10(P.Value),
        fill = logFC,
        color = regulation,
        shape = experiment
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
      scale_shape_manual(values = comparisons, drop = FALSE) +
      scale_colour_manual(values = regulation, drop = FALSE) +
      xlab("Resolvins") + ylab("-log10(P-Value)") +
      guides(color = guide_legend(title = "Comparison")) +
      geom_hline(yintercept = 0 - log10(0.05), linetype = "dashed") +
      guides(
        fill = guide_colourbar(title = "log2(Fold Change)", order = 3),
        color = guide_legend(
          title = "Expression",
          order = 2,
          override.aes = list(shape = 1, size = 4)
        ),
        shape = guide_legend(
          title = "Comparison",
          order = 1,
          override.aes = list(size = 4)
        )
      ) +
      ylim(0, ceiling(max(
        0 - log10(df$P.Value)
      ))) +
      scale_fill_gradientn(
        limits = c(0 - scaleLimits, scaleLimits),
        colours = c(
          "#0000FF",
          "#2E2EFF",
          "#5C5CFF",
          "#8A8AFF",
          "#ffffff",
          "#FF8A8A",
          "#FF5C5C",
          "#FF2E2E",
          "#FF0000"
        )
      )
  )

  print(
    ggplot(
      combinedDf,
      aes(
        x = lipid,
        y = -log10(P.Value),
        fill = logFC,
        color = regulation
      )
    ) +
      geom_point(
        alpha = 1,
        size = 5,
        stroke = 1,
        #position = position_jitter(width = 0.1, height = 0),
        shape = 21
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
      scale_colour_manual(values = regulation, drop = FALSE) +
      xlab("Resolvins") + ylab("-log10(P-Value)") +
      guides(color = guide_legend(title = "Comparison")) +
      geom_hline(yintercept = 0 - log10(0.05), linetype = "dashed") +
      guides(
        fill = guide_colourbar(title = "log2(Fold Change)", order = 3),
        color = guide_legend(
          title = "Expression",
          order = 2,
          override.aes = list(shape = 1, size = 4)
        ),
        shape = guide_legend(
          title = "Comparison",
          order = 1,
          override.aes = list(size = 4)
        )
      ) +
      ylim(0, ceiling(max(
        0 - log10(df$P.Value)
      ))) +
      scale_fill_gradientn(
        limits = c(0 - scaleLimits, scaleLimits),
        colours = c(
          "#0000FF",
          "#2E2EFF",
          "#5C5CFF",
          "#8A8AFF",
          "#ffffff",
          "#FF8A8A",
          "#FF5C5C",
          "#FF2E2E",
          "#FF0000"
        )
      ) +
      facet_wrap(~experiment))
}
