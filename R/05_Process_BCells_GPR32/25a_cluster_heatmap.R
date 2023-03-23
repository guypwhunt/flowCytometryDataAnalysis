try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryName <- "gpr32BCells"
columnNames <- gpr32BCellsClusteringColumnNames

clusterNames <- clusterColumns[3:4]

markersOrCells <- markersOrCellsClassification[3]

markerType <- "Phenotypic"

clusterName <- clusterNames[2]
markersOrCell <- markersOrCells[1]

df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))
df <- as.data.frame(df)
