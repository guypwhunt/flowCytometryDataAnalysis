try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

clusterName <- "clusters_flowsom"

# State if you want to flip the Fold change
flipFoldChange <- FALSE

# Define Signifincance Cut Off
sigCutOff <- 0.05

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)


# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomvisitVisits1, 2ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsombulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsombulbarLimbVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomfastSlowVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomcaseControlVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialStatesStatistics.csv")

##################
recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)


### Differential Abundance ###
DA <- TRUE

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomvisitVisits1, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsomvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomvisitVisits1, 2ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_flowsombulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsombulbarLimbVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomfastSlowVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomcaseControlVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)
