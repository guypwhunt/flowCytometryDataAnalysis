try(source("R/01_functions.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
clusterNames <-
  c(
    "clusters_flowsom",
    "clusters_phenograph",
    "clusters_fast_pg",
    "meta_clusters_flowsom"
  )

markersOrCells <- c("Clusters", "CellPopulations", "Markers")

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    DA <- FALSE


    # State if you want to flip the Fold change
    flipFoldChange <- FALSE

    # Define Signifincance Cut Off
    sigCutOff <- 0.05

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomvisitVisits1, 3AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsomvisitVisits1, 3", markersOrCell, "DifferentialStatesStatistics.csv")
      )

    recalculatePValueAdjustments(DA,
                                 sigCutOff,
                                 fileNames,
                                 clusterName,
                                 markersOrCell,
                                 flipFoldChange)


    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomvisitVisits1, 2AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsomvisitVisits1, 2", markersOrCell, "DifferentialStatesStatistics.csv")
      )


    recalculatePValueAdjustments(DA,
                                 sigCutOff,
                                 fileNames,
                                 clusterName,
                                 markersOrCell,
                                 flipFoldChange)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsombulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsombulbarLimbVisits1", markersOrCell, "DifferentialStatesStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomfastSlowVisits1AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsomfastSlowVisits1", markersOrCell, "DifferentialStatesStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomcaseControlVisits1AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsomcaseControlVisits1", markersOrCell, "DifferentialStatesStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialStatesStatistics.csv",
        paste0("clusters_flowsomcaseControlVisits1, 2, 3", markersOrCell, "DifferentialStatesStatistics.csv")
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)


    ### Differential Abundance ###
    DA <- TRUE

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomvisitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsomvisitVisits1, 3", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )


    recalculatePValueAdjustments(DA,
                                 sigCutOff,
                                 fileNames,
                                 clusterName,
                                 markersOrCell,
                                 flipFoldChange)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomvisitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsomvisitVisits1, 2", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )


    recalculatePValueAdjustments(DA,
                                 sigCutOff,
                                 fileNames,
                                 clusterName,
                                 markersOrCell,
                                 flipFoldChange)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsombulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsombulbarLimbVisits1", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomfastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsomfastSlowVisits1", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomcaseControlVisits1AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsomcaseControlVisits1", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        "clusters_flowsomcaseControlVisits1, 2, 3AllCellsDifferentialAbundanceStatistics.csv",
        paste0("clusters_flowsomcaseControlVisits1, 2, 3", markersOrCell, "DifferentialAbundanceStatistics.csv")
      )

    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)
  }
}
