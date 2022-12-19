try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryNames <- c(#"bCells","monocytes", "senescence", "tCells"
  "gpr18TCells")

directoryName <- "gpr18TCells"

for(directoryName in directoryNames) {
  print(directoryName)
  try({
    df <-
      fread(
        file = paste0(
          "./data/",
          directoryName,
          '/clusteringOutput/phenographClusterStability.csv'
        )#, nrows=100
      )
    df <- as.data.frame(df)
    print(colnames(df))

    df <- df[, !colSums(is.na(df)) == nrow(df)]

    startNumber <- which(colnames(df)=="clusters_phenograph", arr.ind = TRUE)
    endNumber <- ncol(df)

    seqenceNumbner <- seq(from=1, to=endNumber-startNumber)

    colnames(df)[seq(from=startNumber + 1, to=endNumber)] <- paste0("clusters_phenograph", seqenceNumbner)

    fwrite(df, paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/phenographClusterStability.csv'
    ), row.names = FALSE)
  })
}
