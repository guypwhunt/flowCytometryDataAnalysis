try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns

markerName <- "gpr18"

markersOrCells <- markersOrCellsClassification
figureNames <-
  c("DifferentialStatesStatisticscsv",
    "DifferentialAbundanceStatisticscsv")

markersOrCells <- markersOrCells[2]
clusterNames <- clusterNames[3:4]
#figureName <- figureNames[2]

directoryNames <- c("gpr18TCells")

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    for (figureName in figureNames) {
      message("")
      message(clusterName)
      message(markersOrCell)
      message(figureName)
      pattern <-
        paste0(figureName,
               markersOrCell,
               ".csv")
      print(differentialCombinedManhattanPlot(pattern, clusterName, figureName, markersOrCell, directoryNames, markerName))
      #try(dev.off())
    }
  }
}
