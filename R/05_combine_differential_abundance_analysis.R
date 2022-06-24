try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

# State if you want to flip the Fold change
flipFoldChange <- FALSE

# Define Signifincance Cut Off
sigCutOff <- 0.05

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

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
  'Double-Negative T Cells',
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
  'Naive CD8+ T Cells',
  'Naive CD8+ T Cells',
  'Naive CD8+ T Cells',
  'Naive CD8+ T Cells',
  'Naive CD8+ T Cells',
  'Non-Viral Associated Senescent CD8+ T Cells',
  'Non-Viral Associated Senescent CD8+ T Cells',
  'Viral Associated Senescent CD8+ T Cells',
  'Naive CD8+ T Cells',
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

recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames, flipFoldChange)


# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomvisitVisits1, 2ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsombulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsombulbarLimbVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomfastSlowVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomcaseControlVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)


### Differential Abundance ###
DA <- TRUE

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomvisitVisits1, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomvisitVisits1, 2ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsombulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsombulbarLimbVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomfastSlowVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomcaseControlVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, bCellClusterNames, monocyteClusterNames, tCellClusterNames, senescentClusterNames)
