# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
#BiocManager::install("HDCytoData")
#BiocManager::install("CATALYST")

## Load the example dataset
suppressPackageStartupMessages(library(HDCytoData))

# Download and load 'Bodenmiller_BCR_XL' dataset in 'flowSet' format
d_flowSet <- Bodenmiller_BCR_XL_flowSet()

suppressPackageStartupMessages(library(flowCore))

# check data format
d_flowSet

# sample names
pData(d_flowSet)

# dimensions
dim(exprs(d_flowSet[[1]]))

# expression values
exprs(d_flowSet[[1]])[1:6, 1:5]


## Set up meta-data
# Meta-data: experiment information
# check sample order
filenames <- as.character(pData(d_flowSet)$name)

# sample information
sample_id <- gsub("^PBMC8_30min_", "", gsub("//.fcs$", "", filenames))
group_id <- factor(
  gsub("^patient[0-9]+_", "", sample_id), levels = c("Reference", "BCR-XL")
)
patient_id <- factor(gsub("_.*$", "", sample_id))

experiment_info <- data.frame(
  group_id, patient_id, sample_id, stringsAsFactors = FALSE
)
experiment_info

# Meta-data: marker information

# source: Bruggner et al. (2014), Table 1

# column indices of all markers, lineage markers, and functional markers
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# channel and marker names
channel_name <- colnames(d_flowSet)
marker_name <- gsub("//.*$", "", channel_name)

# marker classes
# note: using lineage markers for 'cell type', and functional markers for
# 'cell state'
marker_class <- rep("none", ncol(d_flowSet[[1]]))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

marker_info <- data.frame(
  channel_name, marker_name, marker_class, stringsAsFactors = FALSE
)
marker_info

## Set up design matrix (or model formula)
suppressPackageStartupMessages(library(diffcyt))

# Create design matrix
# note: selecting columns containing group IDs and patient IDs (for an
# unpaired dataset, only group IDs would be included)
design <- createDesignMatrix(
  experiment_info, cols_design = c("group_id", "patient_id")
)

## Set up contrast matrix
# Create contrast matrix
contrast <- createContrast(c(0, 1, rep(0, 7)))

# check
nrow(contrast) == ncol(design)

data.frame(parameters = colnames(design), contrast)


## Differential testing (Option 3: Individual functions)
# Prepare data
d_se <- prepareData(d_flowSet, experiment_info, marker_info)

# Transform data
d_se <- transformData(d_se)

# Generate clusters
# note: include random seed for reproducible clustering
d_se <- generateClusters(d_se, seed_clustering = 123)

# Calculate cluster cell counts
d_counts <- calcCounts(d_se)

# Calculate cluster medians
d_medians <- calcMedians(d_se)

# Test for differential abundance (DA) of cell populations
# Test for differential abundance (DA) of clusters
res_DA <- testDA_edgeR(d_counts, design, contrast)

# display table of results for top DA clusters
topTable(res_DA, format_vals = TRUE)

# calculate number of significant detected DA clusters at 10% false discovery
# rate (FDR)
threshold <- 0.1
table(topTable(res_DA, all = TRUE)$p_adj <= threshold)

# Test for differential states (DS) within cell populations
## UPDATE WITH markers_to_test PARAMETER
# Test for differential states (DS) within clusters
res_DS <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)

# display table of results for top DS cluster-marker combinations
topTable(res_DS, format_vals = TRUE)

# calculate number of significant detected DS cluster-marker combinations at
# 10% false discovery rate (FDR)
threshold <- 0.1
table(topTable(res_DS, all = TRUE)$p_adj <= threshold)


## Differential testing (Option 1)
# Test for differential abundance (DA) of clusters

# note: using default method 'diffcyt-DA-edgeR' and default parameters
# note: include random seed for reproducible clustering
out_DA <- diffcyt(
  d_input = d_flowSet,
  experiment_info = experiment_info,
  marker_info = marker_info,
  design = design,
  contrast = contrast,
  analysis_type = "DA",
  seed_clustering = 123
)

# display table of results for top DA clusters
topTable(out_DA, format_vals = TRUE)

# calculate number of significant detected DA clusters at 10% false discovery
# rate (FDR)
threshold <- 0.1
res_DA_all <- topTable(out_DA, all = TRUE)
table(res_DA_all$p_adj <= threshold)

# Test for differential states (DS) within clusters

# note: using default method 'diffcyt-DS-limma' and default parameters
# note: include random seed for reproducible clustering
out_DS <- diffcyt(
  d_input = d_flowSet,
  experiment_info = experiment_info,
  marker_info = marker_info,
  design = design,
  contrast = contrast,
  analysis_type = "DS",
  seed_clustering = 123,
  plot = FALSE
)

# display table of results for top DS cluster-marker combinations
topTable(out_DS, format_vals = TRUE)

# calculate number of significant detected DS cluster-marker combinations at
# 10% false discovery rate (FDR)
threshold <- 0.1
res_DS_all <- topTable(out_DS, all = TRUE)
table(res_DS_all$p_adj <= threshold)


## Exporting data
# Output object from 'diffcyt()' wrapper function
names(out_DA)
dim(out_DA$d_se)
rowData(out_DA$d_se)
str(assay(out_DA$d_se))
head(assay(out_DA$d_se), 2)

# Extract cell-level table for export as .fcs file
# note: including group IDs, patient IDs, sample IDs, and cluster labels for
# each cell
# note: table must be a numeric matrix (to save as .fcs file)
d_fcs <- assay(out_DA$d_se)
class(d_fcs)

# Save as .fcs file
filename_fcs <- "exported_data.fcs"
write.FCS(
  flowFrame(d_fcs), filename = filename_fcs
)

# Alternatively, save as tab-delimited .txt file
filename_txt <- "exported_data.txt"
write.table(
  d_fcs, file = filename_txt, quote = FALSE, sep = "/t", row.names = FALSE
)


## Visualizations using ‘CATALYST’ package
# Heatmap for top detected DA clusters

# note: use optional argument 'sample_order' to group samples by condition
sample_order <- c(seq(2, 16, by = 2), seq(1, 16, by = 2))

plotHeatmap(out_DA, analysis_type = "DA", sample_order = sample_order)

# Heatmap for top detected DA clusters (alternative code using 'CATALYST')
suppressPackageStartupMessages(library(CATALYST))

# Heatmap for top detected DS cluster-marker combinations

# note: use optional argument 'sample_order' to group samples by condition
sample_order <- c(seq(2, 16, by = 2), seq(1, 16, by = 2))

plotHeatmap(out_DS, analysis_type = "DS", sample_order = sample_order)


###################

df <- read.csv("data/bCells/clusteringOutput/diffusionMapDf.csv")
df <- df[seq_len(10000),]
