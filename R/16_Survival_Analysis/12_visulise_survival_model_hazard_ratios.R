library(dplyr)
library(ggplot2)

figureDirectories <- c("SurvivalAnalysis", "robustSurvivalAnalysis",
                       "ScaledSurvivalAnalysis", "robustScaledSurvivalAnalysis")

for (figureDirectory in figureDirectories) {
  print(figureDirectory)

  figureDirectory <- paste0("data/", figureDirectory, "/cleanResults/")

  fileNames <- list.files(figureDirectory, pattern ="GPR")

  fileNames <- fileNames[!grepl("Univariate", fileNames)]

  fileNames <- append(fileNames, "DiseaseDurationClinicalModel.csv")

  fileNamesPath <- paste0(figureDirectory, fileNames)

  dfs <- lapply(fileNamesPath, read.csv)

  names(dfs) <- fileNames

  i <- 1

  for (df in dfs) {
    df[, "ModelName"] <- names(dfs)[i]

    dfs[[i]] <- df

    i <- i + 1
  }

  df <- do.call(rbind, dfs)

  df <- as.data.frame(df)

  unique(df$ModelName)

  df[df$ModelName == "DiseaseDurationClinicalModel.csv" , "ModelName"] <- "Clinical"
  df[df$ModelName == "GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR32 (FlowSOM)"
  df[df$ModelName == "GPR32clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR32 (Phenograph)"
  df[df$ModelName == "GPR18meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 (FlowSOM)"
  df[df$ModelName == "GPR18clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 (Phenograph)"
  df[df$ModelName == "GPR18GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 & GPR32 (FlowSOM)"
  df[df$ModelName == "GPR18GPR32clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 & GPR32 (Phenograph)"

  df$ModelName <- factor(df$ModelName, levels = c("Clinical", "GPR18 (FlowSOM)", "GPR18 (Phenograph)", "GPR32 (FlowSOM)", "GPR32 (Phenograph)", "GPR18 & GPR32 (FlowSOM)", "GPR18 & GPR32 (Phenograph)"))

  colnames(df)[6] <- "p.value"

  df$lower <- df$Hazard.Ratio - df$Confidence.Interval
  df$upper <- df$Hazard.Ratio + df$Confidence.Interval


  message(paste0("Min: ", min(df$Hazard.Ratio)))
  message(paste0("Max: ", max(df$Hazard.Ratio)))

  dfPositives <- df
  # dfPositives$lower <- dfPositives$Regression.Coefficient - dfPositives$Confidence.Interval
  # dfPositives$upper <- dfPositives$Regression.Coefficient + dfPositives$Confidence.Interval

  dfPositives <- dfPositives[dfPositives$ModelName %in% c("Clinical", "GPR32 (Phenograph)", "GPR18 (Phenograph)", "GPR18 & GPR32 (Phenograph)"), ]

  bCells <- c(
    "F-B (Median GPR18)",
    "F-B (Median GPR32)",
    "CD24+US-B (Median GPR32)",
    "CD24-US-B (Median GPR18)"
  )

  monocytes <- c(
    "HLA-DR-IM-M (Median GPR32)",
    "CD11b Low HLA-DR-aCD11b+IM-M (Median GPR32)",
    "HLA-DR-C-M (Median GPR32)",
    "CD11b Low HLA-DR-aCD11b+C-M (Median GPR32)"
  )

  cells <- list("less", "greater")

  for (cell in cells) {

    if (cell == "less") {
      dfPositivesSub <- dfPositives[dfPositives$Hazard.Ratio < 1,]
    } else {
      dfPositivesSub <- dfPositives[dfPositives$Hazard.Ratio > 1,]
    }

    dfPositivesSub <- dfPositivesSub[dfPositivesSub$Variable %in% c(bCells, monocytes), ]

    dfPositivesSub$Variable <-
      factor(dfPositivesSub$Variable, levels = c(bCells, monocytes))

    print(
      ggplot(
        data = dfPositivesSub,
        aes(
          y = Variable,
          x = Hazard.Ratio,
          xmin = lower,
          xmax = upper
        )
      ) +
        geom_point(size = 3) +
        geom_errorbarh(height = .1) +
        facet_wrap( ~ ModelName) +
        ylab("Cell Populations") +
        xlab("Regression Coefficient [95% Confidence Interval]") +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) #+
        #geom_vline(xintercept = 0, linetype="dotted") #+
        #scale_x_continuous(trans = "log", breaks = scales::breaks_extended())
    )
  }
}

