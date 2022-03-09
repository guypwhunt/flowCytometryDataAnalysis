loadlibraries <- function() {
  library(flowCore)
  library(Biobase)
  library(dplyr)
  library(flowVS)
  library(flowStats)
  library(R.utils)
  library(flowCore)
  library(FlowSOM)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(reshape2)
  library(RColorBrewer)
  library(destiny)
  library(uwot)
  library(slingshot)
  library(cytofkit2)
  library(ggrepel)
  library(tidyverse)
}

preprocessing <- function(directoryName, columnNames, test) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  if (test) {
    filenames <- filenames[1:2]
  }

  ## Defining a function to read a flow cytrometry file in csv format:
  # Each row is a cell, each column is a parameter. In our experience, the flow cytometers sometimes output duplicate entries (listing the same cell twice), we remove these and report.
  # Please check how your csv file is separated and adjust the sep argument in the function if necessary. In this example we import a semicolon separated file.
  read.flow_csv <- function(pathIN){
    raw <- read.csv(pathIN, sep=",", header=TRUE, stringsAsFactors=FALSE)
    IND <- which(duplicated(raw))
    # Check for duplicates and report if found:
    if(any(duplicated(raw))){
      cat(paste0("=== Duplicate entries removed in [",pathIN,"]: ",length(IND)," ===\n"))
      print(head(raw[IND,]))
      cat("----\n")
    }
    return(unique(raw))
  }

  # Read all:
  dfs <- sapply(filenames,read.flow_csv,simplify=FALSE)
  gc()

  ##############################
  #REWRITE TO FLOWFRAME/FLOWSET#
  ##############################

  ## Defining a function to rewrite a csv into a flowframe:
  csv_2_ff <- function(dat){
    # Compute required metadata - column names with description - ranges, min, and max settings
    meta <- data.frame(name=dimnames(dat)[[2]],
                       desc=paste(dimnames(dat)[[2]]),
                       range =(apply(apply(dat,2,range),2,diff)),
                       minRange = apply(dat,2,min),
                       maxRange = apply(dat,2,max))
    # Create flowframe
    flowframef <- new("flowFrame",exprs=as.matrix(dat),parameters=AnnotatedDataFrame(meta))
    return(flowframef)
  }

  # rewrite to flowframe
  dfs_ff = sapply(dfs,function(x) csv_2_ff(x),simplify=FALSE)
  gc()
  rm(dfs)

  # rewrite to flowset
  dfs_fs <- as(dfs_ff,"flowSet")
  gc()
  rm(dfs_ff)

  ###############################
  ####### TRANSFORMATION ########
  ###############################

  ## Each parameter of interest needs to be arcsinh transformed with an individual cofactor. The cofactor can be deduced from the size of the linear region around zero on a biexponential scale, as plotted in a histogram (in conventional gating software).
  # Choose manual transformation or automated transformation (we prefer manual)
  # Define parameters and cofactors for transformations:
  gc()
  jpeg(file = "automatedcofactors.jpeg")
  automatedcofactors <- estParamFlowVS(dfs_fs, columnNames) #this may take a while.
  dev.off()
  gc()

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels=columnNames, cofactor=automatedcofactors)
  gc()
  rm(dfs_fs)
  rm(automatedcofactors)


  ##############################
  ######## NORMALIZATION #######
  ##############################

  ## To correct for technical inter-sample variation we apply normalization by fdaNorm (which automatically detects the number of peaks)
  # We continue with the manual transformed dataset
  # Select the markers which require normalization (based on the densityplots you generated above). Be aware that you don't remove biological variation!
  gc()
  dfs_fs_t_auto_normfda <- warpSet(dfs_fs_t_auto, stains=columnNames)
  gc()

  ##############################
  ####### EXPORT TO FCS ########
  ##############################

  ## The flowset (dfs_fs_t_auto_normfda) can be exported to individual fcs files
  # Create an 'output' folder
  gc()
  dir.create("dataPPOutput", showWarnings = FALSE)
  gc()

  #Save flowframes wihtin flowset as fcs files using the flowCore package
  write.flowSet(dfs_fs_t_auto_normfda, outdir='dataPPOutput', filename = paste0(gsub(".csv", ".fcs", sampleNames(dfs_fs_t_auto_normfda))))
  gc()

  ##############################
  ########### Plots ############
  ##############################

  # Pre-Normalized Plots
  flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

  dir.create("figures", showWarnings = FALSE)
  figureDirectory <- paste0(getwd(),"/figures/")

  if (directoryName == "bCells") {
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotFSCA.jpeg"))
    plot <- densityplot(~FSC.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotSSCA.jpeg"))
    plot <- densityplot(~SSC.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotGPR32AF488A.jpeg"))
    plot <- densityplot(~GPR32...AF488.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotFoxP3PEA.jpeg"))
    plot <- densityplot(~FoxP3.PE.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD19PECF595A.jpeg"))
    plot <- densityplot(~CD19...PE.CF595.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotIgDPerCPCy55A.jpeg"))
    plot <- densityplot(~IgD...PerCP.Cy5.5.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotIFNgPECy7A.jpeg"))
    plot <- densityplot(~IFNg.PE.Cy7.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotFPRL1AF647A.jpeg"))
    plot <- densityplot(~FPRL1...AF647.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotIL.17BV421A.jpeg"))
    plot <- densityplot(~IL.17...BV421.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotZombieNIRA.jpeg"))
    plot <- densityplot(~Zombie.NIR.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD24BV605A.jpeg"))
    plot <- densityplot(~CD24...BV605.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD27BV650.jpeg"))
    plot <- densityplot(~CD27...BV650.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()

    # Normalized Plots
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFSCA.jpeg"))
    plot <- densityplot(~FSC.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotSSCA.jpeg"))
    plot <- densityplot(~SSC.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotGPR32AF488A.jpeg"))
    plot <- densityplot(~GPR32...AF488.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFoxP3PEA.jpeg"))
    plot <- densityplot(~FoxP3.PE.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD19PECF595A.jpeg"))
    plot <- densityplot(~CD19...PE.CF595.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotIgDPerCPCy55A.jpeg"))
    plot <- densityplot(~IgD...PerCP.Cy5.5.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotIFNgPECy7A.jpeg"))
    plot <- densityplot(~IFNg.PE.Cy7.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFPRL1AF647A.jpeg"))
    plot <- densityplot(~FPRL1...AF647.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotIL.17BV421A.jpeg"))
    plot <- densityplot(~IL.17...BV421.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotZombieNIRA.jpeg"))
    plot <- densityplot(~Zombie.NIR.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD24BV605A.jpeg"))
    plot <- densityplot(~CD24...BV605.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD27BV650.jpeg"))
    plot <- densityplot(~CD27...BV650.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
  } else if (directoryName == "senescence") {
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotGPR32.jpeg"))
    plot <- densityplot(~GPR32.AF488.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotKLRG1.jpeg"))
    plot <- densityplot(~KLRG1.PE.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD4.jpeg"))
    plot <- densityplot(~CD4.PE.CF594.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD8PerCP.jpeg"))
    plot <- densityplot(~CD8.PerCP.Cy5.5.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCCR7.jpeg"))
    plot <- densityplot(~CCR7.PE.Cy7.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotFPRL1.jpeg"))
    plot <- densityplot(~FPRL1.AF647.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotZombie.jpeg"))
    plot <- densityplot(~Zombie.NIR.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD27.jpeg"))
    plot <- densityplot(~CD27.BV421.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotx.jpeg"))
    plot <- densityplot(~x.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotPacific.jpeg"))
    plot <- densityplot(~Pacific.Orange.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD45RA.jpeg"))
    plot <- densityplot(~CD45RA.BV605.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD8BV650.jpeg"))
    plot <- densityplot(~CD8.BV650.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"densityPlotCD28.jpeg"))
    plot <- densityplot(~CD28.BV785.A, dfs_fs_t_auto, main="auto")
    print(plot)
    dev.off()
    gc()

    # Normalized Plots
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotGPR32.jpeg"))
    plot <- densityplot(~GPR32.AF488.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotKLRG1.jpeg"))
    plot <- densityplot(~KLRG1.PE.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD4.jpeg"))
    plot <- densityplot(~CD4.PE.CF594.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD8PerCP.jpeg"))
    plot <- densityplot(~CD8.PerCP.Cy5.5.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCCR7.jpeg"))
    plot <- densityplot(~CCR7.PE.Cy7.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFPRL1.jpeg"))
    plot <- densityplot(~FPRL1.AF647.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotZombie.jpeg"))
    plot <- densityplot(~Zombie.NIR.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD27.jpeg"))
    plot <- densityplot(~CD27.BV421.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotx.jpeg"))
    plot <- densityplot(~x.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotPacific.jpeg"))
    plot <- densityplot(~Pacific.Orange.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD45RA.jpeg"))
    plot <- densityplot(~CD45RA.BV605.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD8BV650.jpeg"))
    plot <- densityplot(~CD8.BV650.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD28.jpeg"))
    plot <- densityplot(~CD28.BV785.A, dfs_fs_t_auto_normfda, main="auto")
    print(plot)
    dev.off()
    gc()
  }

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
    setwd("..")
    setwd("..")
  })
  }

convertToDataFrame <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  clinicalData <- read.csv('data/metadata/metadata.csv')

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  ## Optional: when loading clustered fcs files from cytosplore, provide the directory of the text file 'CSPLR_ST.txt'. Cytosplore exports this file upon running the HSNE. This file contains the decoding of the sample numbers.
  pathST <- "X:/Users/guypw/OneDrive/Documents/flowCytometryDataAnalysis/HSNE_clusters_CD4/CSPLR_ST.txt"

  ## Defining a function to read multiple fcs files from a directory 'dir' into a single data.frame:
  # NB: The column in the output named 'fileName' tracks the original file where each cell came from.
  # Optionally perform remapping of column 'CSPLR_ST' holding cytosplore sample numbers to actual names:
  read.flowdat <- function(dir,path_CSPLR_ST=""){
    # Read:
    filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
    flowset <- read.flowSet(files=filepaths, transformation=FALSE, truncate_max_range = FALSE)
    # Transform to data frame:
    x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
    # Map column 'Original' to filename (in this case holding clusters of HSNE):
    filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs", full.names=FALSE))
    names(filenames) <- sort(unique(x$Original))
    x$fileName <- filenames[as.character(x$Original)]
    # Remove column 'Original':
    x <- x[,-which(colnames(x)=="Original")]
    # Optionally remap Cytosplore sample tags to original filename:
    if(file.exists(path_CSPLR_ST)){
      # Read:
      sampID <- gsub(".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[1])))
      names(sampID) <- sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
      x$sampleID <- sampID[as.character(x$CSPLR_ST)]
    }
    return(x)
  }

  ## Read fcs files
  # In our example we will read the data which were clustered in Cytosplore (each fcs file is 1 cluster)
  df <- read.flowdat(dir=dirFCS[1],path_CSPLR_ST = pathST)
  gc()
  write.csv(df, 'dataPPOutput/rawDf.csv')
  gc()

  df <- df[,columnNames]
  gc()
  write.csv(df, 'dataPPOutput/columnsOfInterestDf.csv')
  gc()

  df <- merge(df, clinicalData, by.x = "fileName", by.y = "ï..id")
  df["caseControl"][df["caseControl"] == "Case"] <- 1
  df["caseControl"][df["caseControl"] == "Control"] <- 0

  df["fastSlow"][df["fastSlow"] == "Fast"] <- 1
  df["fastSlow"][df["fastSlow"] == "Slow"] <- 0
  df["fastSlow"][df["fastSlow"] == "N/A"] <- -1
  gc()
  gc()
  write.csv(df, 'dataPPOutput/columnsOfInterestPlusClinicalDataDf.csv')
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}


multipleRegressionTesting <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <- read.csv('dataPPOutput/columnsOfInterestPlusClinicalDataDf.csv')

  # This returns the formula:
  caseModelFormula <- as.formula(paste("caseControl", paste(columnNames, collapse=" + "), sep=" ~ "))

  progressionModelFormula <- as.formula(paste("fastSlow", paste(columnNames, collapse=" + "), sep=" ~ "))

  caseModel <- lm(caseModelFormula, data = df)

  print(summary(caseModel))

  caseModelCoefficient <-
    as.data.frame(summary(caseModel)$coefficient)
  caseModelCoefficient[order(caseModelCoefficient["Pr(>|t|)"]),]

  print(confint(caseModel))

  fastModel <- lm(progressionModelFormula, data = df)

  print(summary(fastModel))

  fastModelCoefficient <-
    as.data.frame(summary(fastModel)$coefficient)
  fastModelCoefficient[order(fastModelCoefficient["Pr(>|t|)"]),]

  print(confint(fastModel))

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}
