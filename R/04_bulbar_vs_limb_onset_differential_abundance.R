bulbarVsLimbOnsetFileNames <-
  c(
    "bulbarVsLimbOnsetOneDifferentialAbundanceStatistics.csv",
    "bulbarVsLimbOnsetVisitOneOneClusterDifferentialAbundanceStatistics.csv",
    "bulbarVsLimbOnsetOneDifferentialStatesStatistics.csv",
    "bulbarVsLimbOnsetVisitOneOneClusterDifferentialStatesStatistics.csv"
  )

directories <- c("bCells", "monocytes")

### B Cells ####
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

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allBCellsDS <-
  read.csv(
    "data/bCells/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialAbundanceStatistics.csv"
  )
allBCellsDS[, "panel"] <-
  "bCells"
allBCellsDS[, "typeOfCells"] <-
  "B Cells"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
bCellsClustersDS <-
  read.csv(
    "data/bCells/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialAbundanceStatistics.csv"
  )
bCellsClustersDS[, "panel"] <-
  "bCells"
bCellsClustersDS[, "typeOfCells"] <- bCellClusterNames

# Combine Datasets
bCellsDF <-
  rbind(
    allBCellsDS,
    bCellsClustersDS
  )

### Monocytes ####
monocyteClusterNames <- c(
  "Non-Classical Monocytes",
  "Classical Monocytes",
  "CD11b+ Classical Monocytes",
  "CD11b+ Intermediate Monocytes",
  "CD11b+ Non-Classical Monocytes",
  "Intermediate Monocytes"
)

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allMonocytesDS <-
  read.csv(
    "data/monocytes/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialAbundanceStatistics.csv"
  )
allMonocytesDS[, "panel"] <-
  "monocytes"
allMonocytesDS[, "typeOfCells"] <-
  "Monocytes"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
monocytesClustersDS <-
  read.csv(
    "data/monocytes/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialAbundanceStatistics.csv"
  )
monocytesClustersDS[, "panel"] <-
  "monocytes"
monocytesClustersDS[, "typeOfCells"] <- monocyteClusterNames

### Combine Datasets ####
monocytesDF <-
  rbind(
    allMonocytesDS,
    monocytesClustersDS
  )

### T Cells ####
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

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allTCellsDS <-
  read.csv(
    "data/tCells/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialAbundanceStatistics.csv"
  )
allTCellsDS[, "panel"] <-
  "tCells"
allTCellsDS[, "typeOfCells"] <-
  "T Cells"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
tCellsClustersDS <-
  read.csv(
    "data/tCells/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialAbundanceStatistics.csv"
  )
tCellsClustersDS[, "panel"] <-
  "tCells"
tCellsClustersDS[, "typeOfCells"] <- tCellClusterNames

### Combine Datasets ####
tCellsDF <-
  rbind(
    allTCellsDS,
    tCellsClustersDS
  )


### Senescent T Cells ####
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

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allSenescenceDS <-
  read.csv(
    "data/senescence/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialAbundanceStatistics.csv"
  )
allSenescenceDS[, "panel"] <-
  "senescence"
allSenescenceDS[, "typeOfCells"] <-
  "Senescent T Cells"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
senescenceClustersDS <-
  read.csv(
    "data/senescence/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialAbundanceStatistics.csv"
  )
senescenceClustersDS[, "panel"] <-
  "senescence"
senescenceClustersDS[, "typeOfCells"] <- senescentClusterNames

### Combine Datasets ####
senescenceDF <-
  rbind(
    allSenescenceDS,
    senescenceClustersDS
  )

### Combine all datasets ###
df <- rbind(bCellsDF, monocytesDF)
df <- rbind(df, tCellsDF)
df <- rbind(df, senescenceDF)

# Bonferroni P-Value Adjustment
df[, "bonferroni_adjusted_p_val"] <- p.adjust(df[, "p_val"], method = "bonferroni")
df[, "minus_log_bonferroni_adjusted_p_val"] <- 0-log10(df[, "bonferroni_adjusted_p_val"])

# Benjamini and Hochberg P-Value Adjustment
df[, "fdr_adjusted_p_val"] <- p.adjust(df[, "p_val"], method = "fdr")
df[, "minus_log_fdr_adjusted_p_val"] <- 0-log10(df[, "fdr_adjusted_p_val"])

# Update differential expression column
df$fdr_diff_expressed <- "NO"
df$fdr_diff_expressed[df[,"logFC"] < 0 & df[,"fdr_adjusted_p_val"] < 0.05] <- "DOWN"
df$fdr_diff_expressed[df[,"logFC"] > 0 & df[,"fdr_adjusted_p_val"] < 0.05] <- "UP"

df$bonferroni_diff_expressed <- "NO"
df$bonferroni_diff_expressed[df[,"logFC"] < 0 & df[,"bonferroni_adjusted_p_val"] < 0.05] <- "DOWN"
df$bonferroni_diff_expressed[df[,"logFC"] > 0 & df[,"bonferroni_adjusted_p_val"] < 0.05] <- "UP"

# Update label
df[,"label"] <- df[,"typeOfCells"]
df$label[is.na(df[,"logFC"])] <- NA
df$label[is.na(df[,"logFC"]) | df[,"fdr_adjusted_p_val"] > 0.05] <- NA

# Define Colours
mycolors <- data.frame(DOWN = "blue",
                       UP = "red",
                       NO = "black")

# GPR32 Plot
par(mar = c(1, 1, 1, 1))
p <- ggplot(data = df,
            aes(
              x = logFC,
              y = minus_log_fdr_adjusted_p_val,
              col = fdr_diff_expressed,
              label = label
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
