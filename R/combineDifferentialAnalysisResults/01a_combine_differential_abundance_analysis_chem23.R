try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

# Define if it is Differential States or Differential Abundance Aanlaysis
clusterNames <- clusterColumns[3:4]

markersOrCells <- markersOrCellsClassification[3]

directories <- c("chem23Monocytes")

markerName <- "chem23"

for(clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
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
                                 markersOrCell, directories, markerName,
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
                                 markersOrCell, directories, markerName,
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
                                 markersOrCell, directories, markerName,
                                 flipFoldChange)


    ## Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "BulbarLimbVisits1AllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "BulbarLimbVisits1",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1LimbAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Limb",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1BulbarAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Bulbar",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)


    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbarAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbar",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowLimbAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1SlowLimb",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastBulbarAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1FastBulbar",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastLimbAllCellsDifferentialStatesStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1FastLimb",
          markersOrCell,
          "DifferentialStatesStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)


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
                                 markersOrCell, directories, markerName,
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
                                 markersOrCell, directories, markerName,
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
                                 markersOrCell, directories, markerName,
                                 flipFoldChange)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "BulbarLimbVisits1AllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "BulbarLimbVisits1",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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


    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

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
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1LimbAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Limb",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1BulbarAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1Bulbar",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)


    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbarAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbar",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastBulbarAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1FastBulbar",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbarAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1SlowBulbar",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

    # Define Directories and files
    fileNames <-
      c(
        paste0(
          clusterName,
          "caseControlVisits1FastLimbAllCellsDifferentialAbundanceStatistics.csv"
        ),
        paste0(
          clusterName,
          "caseControlVisits1FastLimb",
          markersOrCell,
          "DifferentialAbundanceStatistics.csv"
        )
      )

    ##################
    recalculatePValueAdjustments(DA, sigCutOff, fileNames, clusterName, markersOrCell, directories, markerName)

  }
}
