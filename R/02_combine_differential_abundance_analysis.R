# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

# Define Directories and files
directories <- c("bCells", "monocytes", "tCells", "senescence")

fileNames <- c("clusters_flowsomvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

# Define Signifincance Cut Off
sigCutOff <- 0.05

names(fileNames) <- c("allCells", "clusters")

if (DA) {
  replicateValue <- 1
} else {
  replicateValue <- 2
}

### Define Cluster Names ###
bCellClusterNames <- c(
  "Unswitched Memory B Cells",
  "Unswitched Memory B cells",
  "Naive B Cells",
  "Follicular B Cells",
  "Switched Memory B Cells",
  "Immature B Cells",
  "Switched Memory B Cells",
  "Late Memory B Cells"
)

monocyteClusterNames <- c(
  "Non-Classical Monocytes",
  "Classical Monocytes",
  "CD11b+ Classical Monocytes",
  "CD11b+ Intermediate Monocytes",
  "CD11b+ Non-Classical Monocytes",
  "Intermediate Monocytes"
)

tCellClusterNames <- c(
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Memory CD8+ T Cells',
  'Unknown T Cells',
  'Naive CD8+ T Cells',
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Naive T Helper Cells',
  'Naive T Helper Cells',
  'Naive T Helper Cells',
  'Naive T Helper Cells',
  'Memory T Regulatory Cells',
  'Naive T Helper Cells',
  'Naive T Helper Cells',
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Naive T Regulatory Cells',
  'Naive CD8+ T Cells',
  'Memory CD8+ T Cells',
  'Memory CD8+ T Cells',
  'Memory T Helper Cells',
  'Memory T Helper Cells',
  'Memory CD8+ T Cells',
  'Naive CD8+ T Cells',
  'Memory CD8+ T Cells',
  'Memory CD8+ T Cells',
  'Memory CD8+ T Cells'
)

senescentClusterNames <- c(
  'Unknown Senescent T Cells',
  'Unknown Senescent T Cells',
  'Unknown Senescent T Cells',
  'Unknown Senescent T Cells',
  'Unknown Senescent T Cells',
  'Non-Viral Associated Senescent CD8+ T Cells',
  'Non-Viral Associated Senescent CD8+ T Cells',
  'Viral Associated Senescent CD8+ T Cells',
  'Unknown Senescent T Cells',
  'Intermediate Senescent 2',
  'Viral Associated Senescent CD8+ T Cells',
  'Early Senescent',
  'Early Senescent',
  'Viral Associated Senescent CD8+ T Cells',
  'Early Senescent',
  'Intermediate Senescent 2',
  'Viral Associated Senescent CD8+ T Cells',
  'Viral Associated Senescent CD8+ T Cells',
  'Intermediate Senescent 2',
  'Intermediate Senescent 2',
  'Late Senescent',
  'Viral Associated Senescent CD8+ T Cells',
  'Late Senescent',
  'Late Senescent'
)

for (directory in directories) {
  i <- 1
  for (file in fileNames) {
    names(file) <- names(fileNames)[i]
    i <- i +1
    filePath <- paste0("data/", directory, "/differentialTestingOutputs/", file)
    df <- read.csv(filePath)
    df[, "panel"] <- directory
    if (names(file) == "allCells") {
      if (directory == "bCells") {
        df[, "typeOfCells"] <- "B Cells"
      } else if (directory == "monocytes") {
        df[, "typeOfCells"] <- "Monocytes"
      } else if (directory == "tCells") {
        df[, "typeOfCells"] <- "T Cells"
      } else if (directory == "senescence") {
        df[, "typeOfCells"] <- "Senescent T Cells"
      }
    } else if (names(file) == "clusters") {
      if (directory == "bCells") {
        df[, "typeOfCells"] <- rep(bCellClusterNames, replicateValue)
      } else if (directory == "monocytes") {
        df[, "typeOfCells"] <- rep(monocyteClusterNames, replicateValue)
      } else if (directory == "tCells") {
        df[, "typeOfCells"] <- rep(tCellClusterNames, replicateValue)
      } else if (directory == "senescence") {
        df[, "typeOfCells"] <- rep(senescentClusterNames, replicateValue)
      }
    }
    if (exists("combinedDf")) {
      combinedDf <- rbind(combinedDf, df)
    } else {
      combinedDf <- df
    }
  }
}

# Bonferroni P-Value Adjustment
combinedDf[, "bonferroni_adjusted_p_val"] <- p.adjust(combinedDf[, "p_val"], method = "bonferroni")
combinedDf[, "minus_log_bonferroni_adjusted_p_val"] <- 0-log10(combinedDf[, "bonferroni_adjusted_p_val"])

# Benjamini and Hochberg P-Value Adjustment
combinedDf[, "fdr_adjusted_p_val"] <- p.adjust(combinedDf[, "p_val"], method = "fdr")
combinedDf[, "minus_log_fdr_adjusted_p_val"] <- 0-log10(combinedDf[, "fdr_adjusted_p_val"])

# Update differential expression column
combinedDf$bonferroni_diff_expressed <- "NO"
combinedDf$bonferroni_diff_expressed[combinedDf[,"logFC"] < 0 & combinedDf[,"bonferroni_adjusted_p_val"] < sigCutOff] <- "DOWN"
combinedDf$bonferroni_diff_expressed[combinedDf[,"logFC"] > 0 & combinedDf[,"bonferroni_adjusted_p_val"] < sigCutOff] <- "UP"

combinedDf$fdr_diff_expressed <- "NO"
combinedDf$fdr_diff_expressed[combinedDf[,"logFC"] < 0 & combinedDf[,"fdr_adjusted_p_val"] < sigCutOff] <- "DOWN"
combinedDf$fdr_diff_expressed[combinedDf[,"logFC"] > 0 & combinedDf[,"fdr_adjusted_p_val"] < sigCutOff] <- "UP"

# Update labels
combinedDf[,"fdr_label"] <- combinedDf[,"typeOfCells"]
combinedDf$fdr_label[is.na(combinedDf[,"logFC"])] <- NA
combinedDf$fdr_label[is.na(combinedDf[,"logFC"]) | combinedDf[,"fdr_adjusted_p_val"] > sigCutOff] <- NA

combinedDf[,"bonferroni_label"] <- combinedDf[,"typeOfCells"]
combinedDf$bonferroni_label[is.na(combinedDf[,"logFC"])] <- NA
combinedDf$bonferroni_label[is.na(combinedDf[,"logFC"]) | combinedDf[,"bonferroni_adjusted_p_val"] > sigCutOff] <- NA


# Define Colours
mycolors <- data.frame(DOWN = "blue",
                       UP = "red",
                       NO = "black")

dir.create("pValueAdjustmentsResults", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
figureDirectory <- paste0(getwd(), "/figures/")

if (DA) {
  jpeg(file = paste0(
    figureDirectory,
    "fdr",
    str_replace_all(str_replace_all(str_replace_all(fileNames, " ", ""), ",", ""), "\\.", "")[1],
    ".jpeg"
  ))
  par(mar = c(1, 1, 1, 1))
  p <- ggplot(data = combinedDf,
              aes(
                x = logFC,
                y = minus_log_fdr_adjusted_p_val,
                col = fdr_diff_expressed,
                label = fdr_label
              )) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    geom_hline(yintercept = -log10(0.01), col = "red") +
    scale_colour_manual(values = mycolors) +
    geom_text_repel() +
    ggtitle("Differential Abundance of Clusters") +
    xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
  print(p)
  dev.off()
  gc()

  jpeg(file = paste0(
    figureDirectory,
    "bonferroni",
    str_replace_all(str_replace_all(str_replace_all(fileNames, " ", ""), ",", ""), "\\.", "")[1],
    ".jpeg"
  ))
  par(mar = c(1, 1, 1, 1))
  p <- ggplot(data = combinedDf,
              aes(
                x = logFC,
                y = minus_log_bonferroni_adjusted_p_val,
                col = bonferroni_diff_expressed,
                label = bonferroni_label
              )) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    geom_hline(yintercept = -log10(0.01), col = "red") +
    scale_colour_manual(values = mycolors) +
    geom_text_repel() +
    ggtitle("Differential Abundance of Clusters") +
    xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
  print(p)
  dev.off()
  gc()
}
