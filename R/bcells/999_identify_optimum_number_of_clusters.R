try(source("R/01_functions.R"))

loadlibraries()

#################
# variables
directoryName <- "bCells"
columnNames <-
  c("IgD...PerCP.Cy5.5.A",
    "CD24...BV605.A",
    "CD27...BV650.A")

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  apclusterOptimalClusters(directoryName, columnNames)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})

tryCatch({
  veganOptimalClusters(directoryName, columnNames)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})

tryCatch({
  mclustOptimalClusters(directoryName, columnNames)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})
