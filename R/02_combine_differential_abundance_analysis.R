try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
clusterNames <- clusterColumns

# clusterName <- clusterNames[4]
# clusterNames <- clusterNames[3]

markersOrCells <- markersOrCellsClassification

# markersOrCell <- markersOrCells[1]
# markersOrCells <- markersOrCells[3]

my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

foreach(clusterName = clusterNames) %:%
  foreach(markersOrCell = markersOrCells) %dopar% {
    try(source("R/01_functions.R"))

    loadlibraries()

    DA <- FALSE

    # State if you want to flip the Fold change
    flipFoldChange <- FALSE

    # Define Signifincance Cut Off
    sigCutOff <- 0.05

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "visitVisits1, 3AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits1, 3",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
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
        paste0(
          clusterName,
          "visitVisits1, 2AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits1, 2",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
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
        paste0(
          clusterName,
          "visitVisits2, 3AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits2, 3",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA,
                                 sigCutOff,
                                 fileNames,
                                 clusterName,
                                 markersOrCell,
                                 flipFoldChange)


    ## Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "bulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "bulbarLimbVisits1",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "fastSlowVisits1AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "fastSlowVisits1",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Slow",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Fast",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)


    ### Differential Abundance ###
    DA <- TRUE

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "visitVisits1, 3AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits1, 3",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
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
        paste0(
          clusterName,
          "visitVisits1, 2AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits1, 2",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
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
        paste0(
          clusterName,
          "visitVisits2, 3AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "visitVisits2, 3",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
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
        paste0(
          clusterName,
          "bulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "bulbarLimbVisits1",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "fastSlowVisits1AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "fastSlowVisits1",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Slow",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Fast",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell)

  }
