# GPR18 B Cells Inputs
gpr18BCellsColumnNames <- c("IgD", "CD24", "CD27", "GPR18")

gpr18BCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

gpr18BCellsCutoff <- c(0.6, 0.6, 0.6)

# # Chem23 B Cells Inputs
chem23BCellsColumnNames <- c("IgD", "CD24", "CD27", "Chem23")

chem23BCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

chem23BCellsCutoff <- c(0.6, 0.6, 0.6)

# GPR32 B Cells Inputs
gpr32BCellsColumnNames <- c("IgD", "CD24", "CD27", "GPR32")

gpr32BCellsClusteringColumnNames <- c("IgD", "CD24", "CD27")

gpr32BCellsCutoff <- c(0.5, 0.5, 0.5)


# GPR18 Monocyte Inputs
gpr18MonocytesColumnNames <- c("CD11b","CD14", "HLA_DR", "CD16", "CD11b_activated", 'Chem23', "GPR18")

gpr18MonocytesClusteringColumnNames <- c("CD11b","CD14","HLA_DR", "CD16", "CD11b_activated")

gpr18MonocytesCutoff <- c(0.6, 0.5, 0.5, 0.3, 0.4)

# Chem23 Monocyte Inputs
chem23MonocytesColumnNames <- c("CD11b","CD14", "HLA_DR", "CD16", "CD11b_activated", 'GPR18', "Chem23")

chem23MonocytesClusteringColumnNames <- c("CD11b","CD14","HLA_DR", "CD16", "CD11b_activated")

chem23MonocytesCutoff <- c(0.6, 0.5, 0.5, 0.3, 0.4)

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

# Chem23 T Cells
chem23TCellsColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                            "CD45RO", "CD4", "Chem23")

chem23TCellsClusteringColumnNames <- c("CD127", "CD8", "CD25", "FoxP3",
                                      "CD45RO", "CD4")

chem23TCellsCutoff <- c(0.55,0.5,0.4,0.5,0.5,0.4)

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

# chem23 Senescence T Cells
chem23SenescenceColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                "CD8", "CCR7", "Chem23")


chem23SenescenceClusteringColumnNames <- c("CD27", "CD45RA","CD28", "KLRG1", "CD4",
                                          "CD8", "CCR7")

chem23SenescenceCutoff <- c(0.45,0.55,0.55,0.5,0.5,0.6,0.5)


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

# Comparisons

# ALS vs Controls
alsAndSubGroupsVsControlsComparisons <- c(
  "Baseline ALS vs NNC" = 21,
  "Baseline A-F vs NNC"  = 22,
  "Baseline A-S vs NNC" = 23,
  "Baseline A-B vs NNC"  = 24,
  "Baseline A-L vs NNC" = 25
)

progressionSubGroupsVsControlsComparisons <- c(
  "Baseline A-FB vs NNC" = 21,
  "Baseline A-FL vs NNC"  = 22,
  "Baseline A-SB vs NNC" = 23,
  "Baseline A-SL vs NNC"  = 24
  )

progressionSubGroupsVsSubGroupsComparisons <- c(
  "Baseline A-F vs Baseline A-S" = 21,
  "Baseline A-B vs Baseline A-L"  = 22
)

londitudinalComparisons <- c(
  "V2 ALS vs Baseline ALS" = 21,
  "V3 ALS vs Baseline ALS"  = 22,
  "V3 ALS vs V2 ALS" = 23
)

bCellsOrder <- c(
  'B',
  'N-B',
  'LM-B',
  'I-B',
  'F-B',
  'S-B (CD24-)',
  'S-B (CD24+)',
  'US-B (CD24-)',
  'US-B (CD24+)'
)

monocytesOrder <- c(
  'M',
  'NC-M',
  'NC-M (HLA-DR-)',
  'CD11b+NC-M',
  'aCD11b+NC-M (CD11b High)',
  'aCD11b+NC-M (CD11b Low)',
  'aCD11b+NC-M (HLA-DR-, CD11b Low)',
  'IM-M',
  'IM-M (HLA-DR-)',
  'CD11b+IM-M',
  'aCD11b+IM-M (CD11b High)',
  'aCD11b+IM-M (CD11b Low)',
  'aCD11b+IM-M (HLA-DR-, CD11b High)',
  'aCD11b+IM-M (HLA-DR-, CD11b Low)',
  'C-M (HLA-DR-)',
  'CD11b+C-M',
  'aCD11b+C-M (CD11b High)',
  'aCD11b+C-M (HLA-DR-, CD11b High)',
  'aCD11b+C-M (HLA-DR-, CD11b Low)'
)


tCellsOrder <- c(
  'T',
  'N-Tc (CD25-, CD127-, FoxP3+)',
  'N-Tc (CD25+, CD127-, FoxP3-)',
  'N-Tc (CD25+, CD127+, FoxP3+)',
  'M-Tc (CD25-, CD127-, FoxP3-)',
  'M-Tc (CD25-, CD127-, FoxP3+)',
  'M-Tc (CD25-, CD127+, FoxP3+)',
  'M-Tc (CD25+, CD127-, FoxP3-)',
  'M-Tc (CD25+, CD127+, FoxP3-)',
  'M-Tc (CD25+, CD127+, FoxP3+)',
  'N-Tcregs',
  'M-Tcregs',
  'N-Th (CD25-, CD127-, FoxP3-)',
  'N-Th (CD25-, CD127-, FoxP3+)',
  'N-Th (CD25-, CD127+, FoxP3+)',
  'N-Th (CD25+, CD127+, FoxP3+)',
  'M-Th (CD25-, CD127-, FoxP3+)',
  'M-Th (CD25-, CD127+, FoxP3-)',
  'M-Th (CD25-, CD127+, FoxP3+)',
  'M-Th (CD25+, CD127+, FoxP3+)',
  'N-Thregs',
  'M-Thregs',
  'N-Tdp (CD25-, CD127-, FoxP3-)',
  'N-Tdp (CD25-, CD127-, FoxP3+)',
  'N-Tdp (CD25+, CD127-, FoxP3-)',
  'N-Tdp (CD25+, CD127+, FoxP3+)',
  'M-Tdp (CD25-, CD127+, FoxP3+)',
  'M-Tdp (CD25+, CD127+, FoxP3-)',
  'M-Tdp (CD25+, CD127+, FoxP3+)',
  'N-Tdpregs',
  'N-Tdn (CD25-, CD127-, FoxP3-)',
  'N-Tdn (CD25-, CD127-, FoxP3+)',
  'M-Tdn (CD25-, CD127-, FoxP3-)',
  'M-Tdn (CD25-, CD127-, FoxP3+)',
  'M-Tdn (CD25-, CD127+, FoxP3-)',
  'M-Tdn (CD25-, CD127+, FoxP3+)',
  'M-Tdn (CD25+, CD127+, FoxP3-)',
  'M-Tdn (CD25+, CD127+, FoxP3+)',
  'M-Tdnregs'
)

senescenceOrder <- c(
  'S',
  'Viral-S Tc (CD28-)',
  'Viral-S Tc (CD28+)',
  'Viral-S Th',
  'IS2-Tc',
  'LS-Tc',
  'IS1-Tc (KLRG1-)',
  'IS1-Tc (KLRG1+)',
  'ES-Tc',
  'N-Tc (CD27+, CD28-, KLRG1-, CCR7+)',
  'N-Tc (CD27+, CD28+, KLRG1-, CCR7-)',
  'N-Tc (CD27+, CD28+, KLRG1-, CCR7+)',
  'N-Tc (CD27+, CD28+, KLRG1+, CCR+-)',
  'N-Tc (CD27+, CD28+, KLRG1+, CCR7+)',
  'M-Tc (CD27+, KLRG1-, CCR7-, CD28+)',
  'CM-Tc',
  'EM-Tc',
  'N-Th (CD27+, CD28-, KLRG1-, CCR7+)',
  'N-Th (CD27+, CD28+, KLRG1-, CCR7+)',
  'N-Th (CD27+, CD28+, KLRG1+, CCR7+)',
  'M-Th (CD27-, CD28-, KLRG1+, CCR7+)',
  'M-Th (CD27-, CD28+, KLRG1-, CCR7+)',
  'M-Th (CD27-, CD28+, KLRG1+, CCR7+)',
  'M-Th (CD27+, CD28+, KLRG1-, CCR7-)',
  'CM-Th  (KLRG1-, CD28+)',
  'CM-Th (KLRG1-, CD28-)',
  'CM-Th (KLRG1+, CD28+)',
  'EM-Th (CD27-, CD28-, KLRG1+)',
  'EM-Th (CD27-, CD28+, KLRG1+)',
  'EM-Th (CD27+, CD28-, KLRG1+)',
  'M-Tdp',
  'CM-Tdp',
  'EM-Tdp'
)
