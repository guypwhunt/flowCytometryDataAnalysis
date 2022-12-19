try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr18TCells"

markersOrCells <- markersOrCellsClassification

clusterNames <- clusterColumns

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
