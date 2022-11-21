try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

clusterNames <- clusterColumns

markersOrCells <- markersOrCellsClassification
figureNames <-
  c("DifferentialStatesStatisticscsv",
    "DifferentialAbundanceStatisticscsv")

#markersOrCell <- markersOrCells[1]
#clusterName <- clusterNames[2]
#figureName <- figureNames[2]

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
      print(differentialCombinedManhattanPlot(pattern, clusterName, figureName, markersOrCell))
      #try(dev.off())
    }
  }
}
