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
  "Unknown",
  "Follicular B Cells",
  "Switched Memory B Cells",
  "Immature B Cells",
  "Switched Memory B Cells",
  "Late Memory B Cells"
)

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allBCellsDS <-
  read.csv(
    "data/bCells/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialStatesStatistics.csv"
  )
allBCellsDS[, "panel"] <-
  "bCells"
allBCellsDS[, "typeOfCells"] <-
  "B Cells"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
bCellsClustersDS <-
  read.csv(
    "data/bCells/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialStatesStatistics.csv"
  )
bCellsClustersDS[, "panel"] <-
  "bCells"
bCellsClustersDS[, "typeOfCells"] <- c(bCellClusterNames,bCellClusterNames)

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
  "CD11b+ Interstitial Macrophages",
  "CD11b+ Non-Classical Monocytes",
  "Interstitial Macrophages"
)

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allMonocytesDS <-
  read.csv(
    "data/monocytes/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialStatesStatistics.csv"
  )
allMonocytesDS[, "panel"] <-
  "monocytes"
allMonocytesDS[, "typeOfCells"] <-
  "Monocytes"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
monocytesClustersDS <-
  read.csv(
    "data/monocytes/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialStatesStatistics.csv"
  )
monocytesClustersDS[, "panel"] <-
  "monocytes"
monocytesClustersDS[, "typeOfCells"] <- c(monocyteClusterNames,monocyteClusterNames)

### Combine Datasets ####
monocytesDF <-
  rbind(
    allMonocytesDS,
    monocytesClustersDS
  )

### T Cells ####
tCellClusterNames <- c(
  "Non-Classical Monocytes",
  "Classical Monocytes",
  "CD11b+ Classical Monocytes",
  "CD11b+ Interstitial Macrophages",
  "CD11b+ Non-Classical Monocytes",
  "Interstitial Macrophages"
)

# Load B Cell Bulbar Vs Limb Onset Differential States For All B Cells
allTCellsDS <-
  read.csv(
    "data/tCells/differentialTestingOutputs/bulbarVsLimbOnsetVisitOneOneClusterDifferentialStatesStatistics.csv"
  )
allTCellsDS[, "panel"] <-
  "tCells"
allTCellsDS[, "typeOfCells"] <-
  "T Cells"

# Load B Cell Bulbar Vs Limb Onset Differential States For Clusters
TCellsClustersDS <-
  read.csv(
    "data/tCells/differentialTestingOutputs/bulbarVsLimbOnsetOneDifferentialStatesStatistics.csv"
  )
tCellsClustersDS[, "panel"] <-
  "TCells"
tCellsClustersDS[, "typeOfCells"] <- c(tCellClusterNames,tCellClusterNames)

### Combine Datasets ####
tCellsDF <-
  rbind(
    allTCellsDS,
    tCellsClustersDS
  )
##############################################

df <- rbind(monocytesBulbarVsLimbOnsetDf, bCellsBulbarVsLimbOnsetDf)


# Bonferroni P-Value Adjustment
df[, "bonferroni_adjusted_p_val"] <- p.adjust(df[, "p_val"], method = "bonferroni")


# Benjamini and Hochberg P-Value Adjustment
df[, "fdr_adjusted_p_val"] <- p.adjust(df[, "p_val"], method = "fdr")
