try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
DA <- FALSE

clusterName <- "clusters_fast_pg"

# State if you want to flip the Fold change
flipFoldChange <- FALSE

# Define Signifincance Cut Off
sigCutOff <- 0.05

# Define Directories and files
fileNames <- c("clusters_fast_pgvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgvisitVisits1, 3ClustersDifferentialStatesStatistics.csv")

recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)


# Define Directories and files
fileNames <- c("clusters_fast_pgvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgvisitVisits1, 2ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_fast_pgbulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgbulbarLimbVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgfastSlowVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgcaseControlVisits1ClustersDifferentialStatesStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
               "clusters_fast_pgcaseControlVisits1, 2, 3ClustersDifferentialStatesStatistics.csv")

##################
recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)


### Differential Abundance ###
DA <- TRUE

# Define Directories and files
fileNames <- c("clusters_fast_pgvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgvisitVisits1, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_fast_pgvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgvisitVisits1, 2ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, flipFoldChange)

# Define Directories and files
fileNames <- c("clusters_fast_pgbulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgbulbarLimbVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgfastSlowVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgcaseControlVisits1ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)

# Define Directories and files
fileNames <- c("clusters_fast_pgcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
               "clusters_fast_pgcaseControlVisits1, 2, 3ClustersDifferentialAbundanceStatistics.csv")


recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName)
