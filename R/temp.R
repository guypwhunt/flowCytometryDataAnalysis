try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

directories <- c("bCells", "tCells", "monocytes", "senescence")

for (directoryName in directories) {

  if(directoryName == "bCells"){
    columnNames <- bCellsColumnNames
  } else if (directoryName == "tCells"){
    columnNames <- tCellsColumnNames
  } else if (directoryName == "monocytes") {
    columnNames <- monocytesColumnNames
  } else {
    columnNames <- senescenceColumnNames
  }

  df <-
    read.csv(paste0(
      "data/",
      directoryName,
      "/clusteringOutput/",
      "clusteringOutputs.csv"
    ))

  df$GPR32_scaled <- df$GPR32

  rawDf <- read.csv(paste0(
    "data/",
    directoryName,"/dataPPOutput/columnsOfInterestDf.csv"))

  for (col in columnNames) {
    ninetyNinthQuantile <-
      quantile(rawDf[, col], probs = c(0.001, 0.999))
    rawDf <-
      rawDf[rawDf[, col] >= min(ninetyNinthQuantile) &
              rawDf[, col] <= max(ninetyNinthQuantile),]
  }

  message(nrow(df)==nrow(rawDf))

  df$GPR32 <- rawDf$GPR32

  write.csv(df, paste0(
    "data/",
    directoryName,
    "/clusteringOutput/",
    "clusteringOutputs.csv"
  ), row.names = FALSE)
}
