function(directoryName, columnNames, knn) {
  #################
  # variables
  directoryName <- "bCells"
  columnNames <- c("fileName","GPR32...AF488.A","IgD...PerCP.Cy5.5.A",
                   "CD24...BV605.A", "CD27...BV650.A","FPRL1...AF647.A")
  clusterName <- "clusters_flowsom"
  samplesContributionToClustersThreshold <- 10
  differentialAbundanceThreshold <- 0.05
  #################

  # Read experiment data
  experimentInfo <- read.csv("data/metadata/metadata.csv")
  experimentInfo <- experimentInfo[order(experimentInfo[,"sample_id"]),]
  experimentInfo[,"group_id"] <- NA

  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("differentialTestingOutputs", showWarnings = FALSE)

  # Read csv
  df <- read.csv("clusteringOutput/fastPGDf.csv")
  df <- df[order(df[,"fileName"]),]

  # Extract only the relevant columns
  minimalDf <- df[,columnNames]


  # split the dataframe into a dataframe for each file
  listOfDfs <- list()
  for (file in unique(minimalDf[,"fileName"])){
    minimalDfExtract <- minimalDf[minimalDf[,"fileName"] == file,]
    minimalDfExtract <- minimalDfExtract[ , !(names(minimalDfExtract) %in% c("fileName"))]
    listOfDfs <- append(listOfDfs, list(minimalDfExtract))
  }


  # Create marker information
  markerColumnNames <- columnNames[columnNames != "fileName"]
  markerInformation <- data.frame(markerColumnNames)
  markerInformation[,"channel_name"] <- markerColumnNames
  markerInformation[,"marker_name"] <- markerColumnNames
  markerInformation[,"marker_class"] <- rep("type",length(markerColumnNames))
  markerInformation <- markerInformation[,2:ncol(markerInformation)]

  # Transform the input into the correct format
  d_se <- prepareData(listOfDfs, experimentInfo, markerInformation)
  rowData(d_se)[,"cluster_id"] <- df[,clusterName]

  # Calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  rowData(d_counts)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))

  # Transform the cluster cell counts into a plotable format
  counts_df <- assay(d_counts)
  percentage_counts_df <- (counts_df/rowSums(counts_df)) * 100
  t_percentage_counts_df <- t(percentage_counts_df)
  t_percentage_counts_df <- as.data.frame(t_percentage_counts_df)
  colnames(t_percentage_counts_df) <- rowData(d_counts)[,"cluster_id"]
  t_percentage_counts_df[,"rownames"] <- row.names(t_percentage_counts_df)

  figureDirectory <- paste0(getwd(),"/figures/")
  for (columnName in colnames(t_percentage_counts_df)) {
    gc()
    jpeg(file = paste0(figureDirectory,paste0("sampleContributionTo", clusterName, columnName ,".jpeg")))
    par(mar=c(1,1,1,1))

    p <- ggplot(data=t_percentage_counts_df, aes_string(x="rownames",y=columnName)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90)) +
      xlab("Patient Sample") + ylab("Percentage") + ggtitle(columnName)
    print(p)
    dev.off()
    gc()

    try(capture.output(row.names(t_percentage_counts_df[which(t_percentage_counts_df[,columnName]>=samplesContributionToClustersThreshold),]),
                       file = paste0("differentialTestingOutputs/sampleContributionTo", clusterName, columnName ,".txt")))

  }

  # Calculate cluster medians
  d_medians <- calcMedians(d_se)
  rowData(d_medians)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))

  # Create design matrix
  # note: selecting columns containing group IDs and patient IDs (for an
  # unpaired dataset, only group IDs would be included)
  # Updates to case vs control
  experimentInfo[,"group_id"] <- factor(experimentInfo[,"caseControl"])
  experimentInfo[,"patient_id"] <- factor(experimentInfo[,"patient_id"])
  experimentInfo[,"sample_id"] <- factor(experimentInfo[,"sample_id"])
  design <- createDesignMatrix(
    experimentInfo, cols_design = c("group_id")
  )

  # Create contrast (the 1 indicates the columns in the design to test)
  contrast <- createContrast(c(0, 1))

  # Check that design matches control
  nrow(contrast) == ncol(design)
  data.frame(parameters = colnames(design), contrast)

  # Test for differential abundance (DA) of clusters
  res_DA <- testDA_edgeR(d_counts, design, contrast)

  # display table of results for top DA clusters
  topTable(res_DA, format_vals = TRUE)

  # calculate number of significant detected DA clusters at 10% false discovery
  # rate (FDR)
  table(topTable(res_DA, all = TRUE)$p_adj <= differentialAbundanceThreshold)

  # Test for differential states (DS) within clusters
  metadata(d_medians)$id_state_markers <- c(TRUE, rep(FALSE,3), TRUE)
  res_DS <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)

  gc()

  write.csv(as.data.frame(rowData(res_DS)), 'differentialTestingOutputs/differentialStatesStatistics.csv', row.names = FALSE)
  try(saveRDS(res_DS, file = "differentialTestingOutputs/differentialStatesStatistics.rds"))
  write.csv(as.data.frame(rowData(res_DA)), 'differentialTestingOutputs/differentialAbundanceStatistics.csv', row.names = FALSE)
  try(saveRDS(res_DA, file = "differentialTestingOutputs/differentialAbundanceStatistics.rds"))

  rm(res_DA)
  rm(res_DS)
  rm(d_medians)
  rm(df)
  rm(minimalDf)
  rm(listOfDfs)
  rm(markerColumnNames)
  rm(d_se)
  rm(d_counts)
  rm(counts_df)
  rm(t_percentage_counts_df)
  rm(experimentInfo)

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}



# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
#############################################
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(diffcyt))
library(ggplot2)






metadata(d_medians)$id_state_markers

# Experiment
################################
rowRanges(d_counts)
assay(d_counts)
colnames(colData(d_counts))
metadata(d_counts)
rownames(rowData(d_counts))
################################
################################
rowRanges(d_se)
assay(d_se)
colnames(colData(d_se))
metadata(d_se)
rowData(d_se)
################################
################################
rowRanges(d_medians)
assay(d_medians)
colnames(colData(d_medians))
metadata(d_medians)
rowData(d_medians)
################################
################################
rowRanges(res_DA)
assay(res_DA)
colData(res_DA)
metadata(res_DA)
# This contains the stats
rowData(res_DA)
################################
################################
rowRanges(res_DS)
assay(res_DS)
colData(res_DS)
metadata(res_DS)
# This contains the stats
rowData(res_DS)
################################


############################################
# Visit 1 Case vs Control
### Updated to just visit 1 ###
visitOneExperimentInfo <- experimentInfo[experimentInfo[,"visit"] == 1,]
visitOneExperimentInfo <- visitOneExperimentInfo[
  order(visitOneExperimentInfo[,"sample_id"]),]

visitOneMinimalDf <- minimalDf[minimalDf[,"fileName"] %in%
                                 visitOneExperimentInfo[,"sample_id"],]
visitOneMinimalDf <- visitOneMinimalDf[
  order(visitOneMinimalDf[,"fileName"]),]

visitOneDf <- df[df[,"fileName"] %in%
                           visitOneExperimentInfo[,"sample_id"],]
visitOneDf <- visitOneDf[
  order(visitOneDf[,"fileName"]),]

# split the dataframe into a dataframe for each file
listOfDfs <- list()
for (file in unique(visitOneMinimalDf[,"fileName"])){
  minimalDfExtract <- visitOneMinimalDf[visitOneMinimalDf[,"fileName"] == file,]
  minimalDfExtract <- minimalDfExtract[ , !(names(minimalDfExtract) %in% c("fileName"))]
  listOfDfs <- append(listOfDfs, list(minimalDfExtract))
}

# Transform the input into the correct format
d_se <- prepareData(listOfDfs, visitOneExperimentInfo, markerInformation)
head(assay(d_se))
rowData(d_se)[,"cluster_id"] <- visitOneDf[,"clusters_flowsom"]

# Calculate cluster cell counts
d_counts <- calcCounts(d_se)
rowData(d_counts)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))

# Calculate cluster medians
d_medians <- calcMedians(d_se)
rowData(d_medians)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))
assay(d_medians)

# Create design matrix
# note: selecting columns containing group IDs and patient IDs (for an
# unpaired dataset, only group IDs would be included)
# Updates to case vs control
visitOneExperimentInfo[,"group_id"] <- factor(visitOneExperimentInfo[,"caseControl"])
visitOneExperimentInfo[,"patient_id"] <- factor(visitOneExperimentInfo[,"patient_id"])
visitOneExperimentInfo[,"sample_id"] <- factor(visitOneExperimentInfo[,"sample_id"])
visitOneExperimentInfo[,"gender"] <- factor(visitOneExperimentInfo[,"gender"])
visitOneExperimentInfo[,"ageAtVisit"] <- factor(visitOneExperimentInfo[,"ageAtVisit"])
visitOneExperimentInfo[,"ageAtVisitDouble"] <- as.double(visitOneExperimentInfo[,"ageAtVisit"])

head(visitOneExperimentInfo)
design <- createDesignMatrix(
  visitOneExperimentInfo, cols_design = c("group_id"
  ,"gender","ageAtVisitDouble")
)

# Create contrast (the 1 indicates the columns in the design to test)
contrast <- createContrast(c(0, 1, 0, 0)
                           #, rep(0, 89)
)

# Check that design matches control
nrow(contrast) == ncol(design)
data.frame(parameters = colnames(design), contrast)

# Test for differential abundance (DA) of clusters
res_DA <- testDA_edgeR(d_counts, design, contrast)

# display table of results for top DA clusters
topTable(res_DA, format_vals = TRUE)

# calculate number of significant detected DA clusters at 10% false discovery
# rate (FDR)
table(topTable(res_DA, all = TRUE)$p_adj <= threshold)

# Test for differential states (DS) within clusters
metadata(d_medians)$id_state_markers <- c(TRUE, rep(FALSE,3), TRUE)
res_DS <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)

# display table of results for top DS cluster-marker combinations
resultsTable <- topTable(res_DS, format_vals = TRUE)
as.data.frame(resultsTable)

# calculate number of significant detected DS cluster-marker combinations at
# 10% false discovery rate (FDR)
table(topTable(res_DS, all = TRUE)$p_adj <= threshold)


############################################
# Visit 1 Fast vs Slow Progression
### Updated to just visit 1 ###
visitOneExperimentInfo <- experimentInfo[experimentInfo[,"visit"] == 1,]
visitOneExperimentInfo <- visitOneExperimentInfo[visitOneExperimentInfo[,"caseControl"] == "Case",]
visitOneExperimentInfo <- visitOneExperimentInfo[
  order(visitOneExperimentInfo[,"sample_id"]),]
head(visitOneExperimentInfo)

visitOneMinimalDf <- minimalDf[minimalDf[,"fileName"] %in%
                                 visitOneExperimentInfo[,"sample_id"],]
visitOneMinimalDf <- visitOneMinimalDf[
  order(visitOneMinimalDf[,"fileName"]),]

visitOneDf <- df[df[,"fileName"] %in%
                   visitOneExperimentInfo[,"sample_id"],]
visitOneDf <- visitOneDf[
  order(visitOneDf[,"fileName"]),]

# split the dataframe into a dataframe for each file
listOfDfs <- list()
for (file in unique(visitOneMinimalDf[,"fileName"])){
  minimalDfExtract <- visitOneMinimalDf[visitOneMinimalDf[,"fileName"] == file,]
  minimalDfExtract <- minimalDfExtract[ , !(names(minimalDfExtract) %in% c("fileName"))]
  listOfDfs <- append(listOfDfs, list(minimalDfExtract))
}

# Transform the input into the correct format
d_se <- prepareData(listOfDfs, visitOneExperimentInfo, markerInformation)
head(assay(d_se))
rowData(d_se)[,"cluster_id"] <- visitOneDf[,"clusters_flowsom"]

# Calculate cluster cell counts
d_counts <- calcCounts(d_se)
rowData(d_counts)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))

# Calculate cluster medians
d_medians <- calcMedians(d_se)
rowData(d_medians)[,"cluster_id"] <- as.factor(rownames(rowData(d_counts)))
assay(d_medians)

# Create design matrix
# note: selecting columns containing group IDs and patient IDs (for an
# unpaired dataset, only group IDs would be included)
# Updates to case vs control
visitOneExperimentInfo[,"group_id"] <- factor(visitOneExperimentInfo[,"fastSlow"])
visitOneExperimentInfo[,"patient_id"] <- factor(visitOneExperimentInfo[,"patient_id"])
visitOneExperimentInfo[,"sample_id"] <- factor(visitOneExperimentInfo[,"sample_id"])
visitOneExperimentInfo[,"gender"] <- factor(visitOneExperimentInfo[,"gender"])
visitOneExperimentInfo[,"ageAtVisit"] <- factor(visitOneExperimentInfo[,"ageAtVisit"])
visitOneExperimentInfo[,"ageAtVisitDouble"] <- as.double(visitOneExperimentInfo[,"ageAtVisit"])
visitOneExperimentInfo[,"bulbarLimb"] <- factor(visitOneExperimentInfo[,"bulbarLimb"])

head(visitOneExperimentInfo)
design <- createDesignMatrix(
  visitOneExperimentInfo, cols_design = c("group_id"
                                          ,"gender","ageAtVisitDouble",
                                          "bulbarLimb")
)

# Create contrast (the 1 indicates the columns in the design to test)
contrast <- createContrast(c(0, 1, 0, 0, 0)
                           #, rep(0, 89)
)

# Check that design matches control
nrow(contrast) == ncol(design)
data.frame(parameters = colnames(design), contrast)

# Test for differential abundance (DA) of clusters
res_DA <- testDA_edgeR(d_counts, design, contrast)

# display table of results for top DA clusters
topTable(res_DA, format_vals = TRUE)

# calculate number of significant detected DA clusters at 10% false discovery
# rate (FDR)
table(topTable(res_DA, all = TRUE)$p_adj <= threshold)

# Test for differential states (DS) within clusters
metadata(d_medians)$id_state_markers <- c(TRUE, rep(FALSE,3), TRUE)
res_DS <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)

# display table of results for top DS cluster-marker combinations
resultsTable <- topTable(res_DS, format_vals = TRUE)
as.data.frame(resultsTable)

# calculate number of significant detected DS cluster-marker combinations at
# 10% false discovery rate (FDR)
table(topTable(res_DS, all = TRUE)$p_adj <= threshold)

