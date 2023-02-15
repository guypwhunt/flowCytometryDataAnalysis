try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns

markerNames <- c("gpr18", "gpr32")

markersOrCells <- markersOrCellsClassification
figureNames <-
  c("DifferentialStatesStatisticscsv",
    "DifferentialAbundanceStatisticscsv"
  )

markersOrCells <- markersOrCells[3]
clusterNames <- clusterNames[3:4]

directoryNames <- c(
  "BCells",
  "Monocytes",
  "TCells",
  "Senescence"
)

# clusterName <- clusterNames[2]
# markersOrCell <- markersOrCells[1]
# figureName <- figureNames[1]

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    for (figureName in figureNames) {
      message("")
      message(clusterName)
      message(markersOrCell)
      message(figureName)
      pattern <-
        paste0(figureName,
               markersOrCell,
               ".csv")

      directory <- paste0("data/", markerNames, "pValueAdjustmentsResults/")

      fileNames <- list.files(directory, pattern = pattern)

      fileNames <- fileNames[grepl(clusterName, fileNames, fixed = TRUE)]

      fileNamesPath <- paste0(directory, fileNames)

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

      for (directoryName in directoryNames) {
        for (markerName in markerNames) {
          try({
            stability <-
              fread(
                paste0(
                  "data/",
                  markerName,
                  directoryName,
                  "/clusteringOutput/",
                  clusterName,
                  "MarkerStability",
                  ".csv"
                )
              )
            stability <- as.data.frame(stability)
            numberOfColumns <- ncol(stability)

            stability$marker <- toupper(markerName)

            stability$median <- NA

            for (row in seq(nrow(stability))) {
              stability[row, "median"] <-
                median(as.double(stability[row, seq(from = numberOfColumns - 99, to = numberOfColumns)]))
            }

            clustermarkerName <- paste0(clusterName, "Marker")

            try({
              if (!exists("combinedStability")) {
                combinedStability <- stability[, c(clustermarkerName, "marker", "median")]
              } else {
                combinedStability <-
                  rbind(combinedStability, stability[, c(clustermarkerName, "marker", "median")])
              }
            })
          })
        }
      }

      combinedDf[grepl("gpr18", combinedDf$panel), "Marker"] <- "GPR18"
      combinedDf[grepl("gpr32", combinedDf$panel), "Marker"] <- "GPR32"
      combinedDf[grepl("chem23", combinedDf$panel), "Marker"] <- "Chem23"


      try(combinedDf <- merge(combinedDf, combinedStability, by.x = c("cluster_id", "Marker"), by.y = c(clustermarkerName, "marker"), all.x = TRUE))

      combinedDf[combinedDf$cluster_id == 1, "median"] <- 1

      combinedDf$jaccard <- NA
      combinedDf$shape <- 18

      # nrow(combinedDf)

      combinedDf <- combinedDf[combinedDf$panel %in% apply(expand.grid(markerNames, directoryNames), 1, paste, collapse=""),]

      combinedDf[is.na(combinedDf$median), "median"] <- 0

      combinedDf[combinedDf$median < 0.65, "jaccard"] <- "Low"
      combinedDf[combinedDf$median > 0.65, "jaccard"] <- "Moderate"
      combinedDf[combinedDf$median > 0.75, "jaccard"] <- "High"
      combinedDf[combinedDf$median > 0.85, "jaccard"] <- "Very High"

      combinedDf[combinedDf$jaccard == "Low", "shape"] <- 18
      combinedDf[combinedDf$jaccard == "Moderate", "shape"] <- 16
      combinedDf[combinedDf$jaccard == "High", "shape"] <- 15
      combinedDf[combinedDf$jaccard == "Very High", "shape"] <- 17

      combinedDf <- combinedDf[order(combinedDf$median),]

      combinedDf$jaccard <- factor(combinedDf$jaccard, levels = c("Very High", "High", "Moderate", "Low"))

      maxNumberOfComparisons <- 0

      combinedDf <- combinedDf[!duplicated(combinedDf[, c(
        "cluster_id",
        "typeOfCells",
        "panel",
        "Marker",
        "experiment",
        "logFC",
        "logCPM",
        "AveExpr",
        "t",
        "B",
        "LR",
        "p_val"
      )]),]

      for (comparison in unique(combinedDf$experiment)) {
        if (figureName == "DifferentialStatesStatisticscsv") {
          x <- combinedDf$experiment == comparison
          combinedDf[x, "fdr_adjusted_p_val"] <-
            p.adjust(combinedDf[x, "p_val"], method = "fdr")
          y <- nrow(combinedDf[x, ])
          if (y > maxNumberOfComparisons) {
            maxNumberOfComparisons <- y
          }
        } else {
          for (marker in unique(combinedDf$Marker)) {
            x <- combinedDf$experiment == comparison & combinedDf$Marker == marker
            combinedDf[x, "fdr_adjusted_p_val"] <-
              p.adjust(combinedDf[x, "p_val"], method = "fdr")
            y <- nrow(combinedDf[x, ])
            if (y > maxNumberOfComparisons) {
              maxNumberOfComparisons <- y
            }
          }
        }
      }

      combinedDf$minus_log_fdr_adjusted_p_val <- 0 - log10(combinedDf$fdr_adjusted_p_val)
      combinedDf$minus_log_p_val <- 0 - log10(combinedDf$p_val)

      combinedDf <- updateMarkerNames(combinedDf)

      combinedDf[combinedDf$experiment == paste0(clusterName,"BulbarLimbVisits1AllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Bulbar vs Limb"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1AllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "ALS vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1BulbarAllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Bulbar vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastAllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Fast vs Controls"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastBulbarAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast Bulbar vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastLimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast Limb vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1LimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Limb vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowBulbarAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow Bulbar vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowLimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow Limb vs Control"
      combinedDf[combinedDf$experiment == paste0(clusterName,"fastSlowVisits1AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast vs Slow"
      combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits12AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 2 vs Visit 1"
      combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits13AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 3 vs Visit 1"
      combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits23AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 3 vs Visit 2"

      combinedDf$panel <- factor(combinedDf$panel)
      combinedDf <- combinedDf[order(combinedDf$typeOfCells),]
      combinedDf <- combinedDf[order(combinedDf$panel),]

      combinedDf$typeOfCells <- factor(combinedDf$typeOfCells, levels =  unique(combinedDf$typeOfCells))

      try({
        combinedDf <- combinedDf[, c(
          "cluster_id",
          "typeOfCells",
          "panel",
          "Marker",
          "experiment",
          "logFC",
          "logCPM",
          "AveExpr",
          "t",
          "B",
          "LR",
          "p_val",
          "fdr_adjusted_p_val",
          "median",
          "minus_log_p_val",
          "minus_log_fdr_adjusted_p_val"
        )]

        colnames(combinedDf) <- c(
          "ID",
          "Cell Population Name",
          "Panel",
          "Marker",
          "Comparison",
          "Log Fold Change",
          "Log Counts Per Million",
          "Average Expression",
          "T Statistic",
          "B Statistic (Log Odds)",
          "Likelihood-Ratio",
          "Raw P-Value",
          "FDR Adjusted P-Value",
          "Jaccard Similarity Coefficient",
          "Minus Log (Raw P-Value)",
          "Minus Log (FDR Adjusted P-Value)"
        )
      })

      combinedDf <- filterFlowCytometryComparisons(combinedDf)

      dir.create("data/pValueAdjustmentsResults", showWarnings = FALSE)

      fwrite(combinedDf,
             paste0("data/pValueAdjustmentsResults/", clusterName, markersOrCell, figureName, ".csv"))

      rm(list=ls()[!ls() %in% c("clusterNames", "markerNames", "markersOrCells",
                     "figureNames", "directoryNames",
                     "clusterName", "markersOrCell", "figureName")])

      try(source("R/01_functions.R"))
      try(source("R/00_datasets.R"))
    }
  }
}
