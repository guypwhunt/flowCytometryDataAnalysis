try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

clusterName <- "clusters_phenograph"

# State if you want to flip the Fold change
flipFoldChange <- FALSE

# Define Signifincance Cut Off
sigCutOff <- 0.05

# Define Directories and files
fileNames <- c("clusters_phenographvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)


# Define Directories and files
fileNames <- c("clusters_phenographvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographvisitVisits1, 2ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_phenographbulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographbulbarLimbVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographfastSlowVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographcaseControlVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_phenographcaseControlVisits1, 2, 3ClustersDifferentialStatesStatistics.csv")

##################
recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)


### Differential Abundance ###
DA <- TRUE

# Define Directories and files
fileNames <- c("clusters_phenographvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographvisitVisits1, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_phenographvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographvisitVisits1, 2ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_phenographbulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographbulbarLimbVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographfastSlowVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographcaseControlVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_phenographcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_phenographcaseControlVisits1, 2, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)
