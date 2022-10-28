# B Cells Inputs
bCellsColumnNames <- c("IgD", "CD24", "CD27", "GPR32")

bCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

# Monocyte Inputs
monocytesColumnNames <- c("CD11b","CD14", "HLA_DR", "CD16", "CD11b_activated", "GPR32")

monocytesClusteringColumnNames <- c("CD11b","CD14","HLA_DR", "CD16", "CD11b_activated")

# T Cells Inputs
tCellsColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                       "CD45RO", "CD4", "GPR32")


tCellsClusteringColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                                 "CD45RO", "CD4")

# Senescence T Cells
senescenceColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                           "CD8", "CCR7", "GPR32")


senescenceClusteringColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                     "CD8", "CCR7")


# Clustering Inputs
numberOfClusters <- seq(3,20)
knn <- 50

# Visulisisation Inputs
clusterColumns <- c(
  #"clusters_flowsom",
  #"clusters_fast_pg",
  "meta_clusters_flowsom",
  "clusters_phenograph"
)
markersOrCellsClassification <- c(#"Clusters",
  #"CellPopulations",
  "Markers")

# Parallel inputs
n.cores <- 2

# Cluster Cutoff
cutOff <- 0.4

iterations <- 100
