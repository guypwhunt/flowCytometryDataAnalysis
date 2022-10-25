try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns
markersOrCells <- markersOrCellsClassification
figureNames <-
  c("DifferentialStatesStatisticscsv",
    "DifferentialAbundanceStatisticscsv")

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    for (figureName in figureNames) {
      pattern <-
        paste0(figureName,
               markersOrCell,
               ".csv")
      differentialCombinedManhattanPlot(pattern, clusterName, figureName, markersOrCell)
      try(dev.off())
    }
  }
}
