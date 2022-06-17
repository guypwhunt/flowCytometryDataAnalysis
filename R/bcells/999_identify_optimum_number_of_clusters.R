try(source("R/01_functions.R"))

loadlibraries()

#################
# variables
directoryName <- "bCells"
columnNames <-
  c("IgD...PerCP.Cy5.5.A",
    "CD24...BV605.A",
    "CD27...BV650.A")


minClusters <- 1
maxClusters <- 3

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  veganOptimalClusters(directoryName, columnNames, minClusters, maxClusters)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})


minClusters <- 3
maxClusters <- 5

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  veganOptimalClusters(directoryName, columnNames, minClusters, maxClusters)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})

minClusters <- 5
maxClusters <- 7

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  veganOptimalClusters(directoryName, columnNames, minClusters, maxClusters)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})

minClusters <- 7
maxClusters <- 9

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  veganOptimalClusters(directoryName, columnNames, minClusters, maxClusters)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})

minClusters <- 9
maxClusters <- 11

# https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
tryCatch({
  veganOptimalClusters(directoryName, columnNames, minClusters, maxClusters)
},
error = function(cond) {
  setwd("..")
  setwd("..")
})
