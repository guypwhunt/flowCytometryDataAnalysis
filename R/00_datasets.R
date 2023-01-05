# GPR18 B Cells Inputs
gpr18BCellsColumnNames <- c("IgD", "CD24", "CD27", "GPR18")

gpr18BCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

gpr18BCellsCutoff <- c(0.6, 0.6, 0.6)

# GPR32 B Cells Inputs
gpr32BCellsColumnNames <- c("IgD", "CD24", "CD27", "GPR32")

gpr32BCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

gpr32BCellsCutoff <- c(0.5, 0.5, 0.5)

# GPR18 Monocyte Inputs
gpr18MonocytesColumnNames <- c("CD11b","CD14", "HLA_DR", "CD16", "CD11b_activated", 'Chem23', "GPR18")

gpr18MonocytesClusteringColumnNames <- c("CD11b","CD14","HLA_DR", "CD16", "CD11b_activated")

gpr18MonocytesCutoff <- c(0.6, 0.5, 0.5, 0.3, 0.4)

# GPR32 Monocyte Inputs
gpr32MonocytesColumnNames <- c("CD11b","CD14", "HLA_DR", "CD16", "CD11b_activated", "GPR32")

gpr32MonocytesClusteringColumnNames <- c("CD11b","CD14","HLA_DR", "CD16", "CD11b_activated")

gpr32MonocytesCutoff <- c(0.6, 0.5, 0.5, 0.35, 0.5)

# GPR18 T Cells
gpr18TCellsColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                            "CD45RO", "CD4", "GPR18")

gpr18TCellsClusteringColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                                      "CD45RO", "CD4")

gpr18TCellsCutoff <- c(0.55,0.5,0.4,0.5,0.5,0.4)

# GPR32 T Cells Inputs
gpr32TCellsColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                       "CD45RO", "CD4", "GPR32")


gpr32TCellsClusteringColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                                 "CD45RO", "CD4")

gpr32TCellsCutoff <- c(0.5,0.6,0.5,0.4,0.5,0.6)

# GPR18 Senescence T Cells
gpr18SenescenceColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                "CD8", "CCR7", "GPR18")


gpr18SenescenceClusteringColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                          "CD8", "CCR7")

gpr18SenescenceCutoff <- c(0.45,0.55,0.55,0.5,0.5,0.6,0.5)

# GPR32 Senescence T Cells
gpr32SenescenceColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                "CD8", "CCR7", "GPR32")


gpr32SenescenceClusteringColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                          "CD8", "CCR7")

gpr32SenescenceCutoff <- c(0.5,0.5,0.6,0.5,0.5,0.6,0.5)

# Clustering Inputs
numberOfClusters <- seq(3,20)
knn <- 50

# Visulisisation Inputs
clusterColumns <- c(
  "clusters_flowsom",
  "clusters_fast_pg",
  "meta_clusters_flowsom",
  "clusters_phenograph"
)
markersOrCellsClassification <- c(
  "Clusters",
  "CellPopulations",
  "Markers"
  )

# Parallel inputs
n.cores <<- 11

# Cluster Cutoff
cutOff <- 0.4

iterations <- 100

clusterStabilityCutoff <- 0.85
