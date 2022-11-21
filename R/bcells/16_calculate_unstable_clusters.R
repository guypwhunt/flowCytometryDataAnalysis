try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

#clusterName <- clusterColumns[1]

directoryName <- "bCells"

markersOrCells <- markersOrCellsClassification

#markersOrCell <- markersOrCells[1]

clusterNames <- clusterColumns

#clusterName <- clusterNames[3]

cutOff <- cutOff

for (clusterName in clusterNames) {
  for (markersOrCell in markersOrCells) {
    df <- fread(file=
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        markersOrCell,
        "Counts.csv"
      )
    )
	
	df <- as.data.frame(df)

    identifyUnstableClustersFromCounts(directoryName, markersOrCell,
                    clusterName, df, cutOff)
  }
}
