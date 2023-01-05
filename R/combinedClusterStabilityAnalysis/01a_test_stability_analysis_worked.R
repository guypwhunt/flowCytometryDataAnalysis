try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryNames <- c(
  # "gpr18BCells",
  # "gpr18Monocytes",
  # "gpr18Senescence",
  # "gpr18TCells",
  # "gpr32BCells",
  # "gpr32Monocytes",
  "gpr32Senescence",
  "gpr32TCells"
)

for (directoryName in directoryNames) {
  try({
    clusteringOutputs <-
      fread(paste0(
        'data/',
        directoryName,
        '/clusteringOutput/clusteringOutputs.csv'
      ))
    clusteringOutputs <- as.data.frame(clusteringOutputs)

    phenographClusterStability <- fread(
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/phenographClusterStability.csv"
      )
    )
    phenographClusterStability <-
      as.data.frame(phenographClusterStability)

    if(nrow(clusteringOutputs) == nrow(phenographClusterStability)){
      fwrite(phenographClusterStability, paste0(
        "data/",
        directoryName,
        "/clusteringOutput/phenographClusterStability-backup.csv"
      ))
    }

    message(directoryName)
    message(nrow(clusteringOutputs) == nrow(phenographClusterStability))
    message(paste0("clusteringOutputs: ", nrow(clusteringOutputs)))
    message(paste0(
      "phenographClusterStability: ",
      nrow(phenographClusterStability)
    ))
    message(paste0(
      "phenographClusterStability ncol: ",
      ncol(phenographClusterStability)
    ))

    message(paste0(
      "non-na phenographClusterStability ncol: ",
      ncol(phenographClusterStability) - sum(colSums(is.na(phenographClusterStability)) == ncol(phenographClusterStability), na.rm = TRUE)
    ))

    # for (n in seq(ncol(clusteringOutputs))){
    #   try({
    #     message("")
    #     message(colnames(clusteringOutputs)[n])
    #     message(paste0("Max: ", max(clusteringOutputs[,n])))
    #     message(paste0("Min: ", min(clusteringOutputs[,n])))
    #   }
    #   )
    # }
    print(colnames(phenographClusterStability)[ncol(phenographClusterStability)])
    message("")
    rm(clusteringOutputs)
    rm(phenographClusterStability)
    gc()
  })
}


# for (directoryName in directoryNames) {
#   try({
#     message(directoryName)
#     clusteringOutputs <-
#       fread(paste0(
#         'data/',
#         directoryName,
#         '/clusteringOutput/clusteringOutputs.csv'
#       ))
#     clusteringOutputs <- as.data.frame(clusteringOutputs)
#     print(unique(clusteringOutputs$fileName))
#     message("")
#   })
# }

# experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")
#
# experimentInfo <- as.data.frame(experimentInfo)
#
# experimentInfo <- experimentInfo[experimentInfo$experiment == "flowCytometry", ]
#
# experimentInfo <- experimentInfo[, c("patient_id", "sample_id")]
#
# for (directoryName in directoryNames) {
#   try({
#     message(directoryName)
#     clinical <- updateClinicalData(experimentInfo, directoryName)
#     clusteringOutputs <-
#       as.data.frame(fread(paste0(
#         'data/',
#         directoryName,
#         '/clusteringOutput/clusteringOutputs.csv'
#       )))
#
#     fileNames <- unique(clusteringOutputs$fileName)
#
#     fileNames <- fileNames[!fileNames %in% clinical$sample_id]
#     print(fileNames)
#
#     message("")
#   })
# }
