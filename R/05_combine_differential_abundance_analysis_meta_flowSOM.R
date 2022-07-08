try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

clusterName <- "meta_clusters_flowsom"

# State if you want to flip the Fold change
flipFoldChange <- FALSE

# Define Signifincance Cut Off
sigCutOff <- 0.05

# Define Directories and files
fileNames <- c("meta_clusters_flowsomvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsomvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)


# Define Directories and files
fileNames <- c("meta_clusters_flowsomvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsomvisitVisits1, 2ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("meta_clusters_flowsombulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsombulbarLimbVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsomfastSlowVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsomcaseControlVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
               "meta_clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialStatesStatistics.csv")

##################
recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)


### Differential Abundance ###
DA <- TRUE

# Define Directories and files
fileNames <- c("meta_clusters_flowsomvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsomvisitVisits1, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsomvisitVisits1, 2ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("meta_clusters_flowsombulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsombulbarLimbVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsomfastSlowVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsomcaseControlVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("meta_clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
               "meta_clusters_flowsomcaseControlVisits1, 2, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)
