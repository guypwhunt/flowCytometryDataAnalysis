loadlibraries <- function() {
  library(R.utils)
  library(stringr)
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
  library(diffcyt)
}

ungzipFiles <- function() {
  workingDirectory <- getwd()

  dataDirectorys <- c("/data/bCells"
                      ,"/data/monocytes",
                      "/data/senescence","/data/tCells"
  )

  for (directory in dataDirectorys) {
    try(setwd(paste0(workingDirectory,directory)))

    gzFilenames <- list.files(pattern = ".gz")

    for (gzFilename in gzFilenames) {
      gunzip(gzFilename)
    }
  }

  setwd(workingDirectory)
  gc()
}

gzipFiles <- function() {
  workingDirectory <- getwd()

  dataDirectorys <- c("/data/bCells","/data/monocytes",
                      "/data/senescence","/data/tCells")

  for (directory in dataDirectorys) {
    try(setwd(paste0(workingDirectory,directory)))

    print(getwd())

    filenames <- list.files(pattern = ".csv")

    for (filename in filenames) {
      gzip(filename)
    }
  }

  setwd(workingDirectory)

  gc()
}

preprocessing <- function(directoryName, columnNames, test) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  # Create an 'output' folder
  gc()
  dir.create("dataPPOutput", showWarnings = FALSE)
  gc()

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  if (test) {
    filenames <- filenames[1:2]
  }

  ## Defining a function to read a flow cytrometry file in csv format:
  # Each row is a cell, each column is a parameter. In our experience, the
  # flow cytometers sometimes output duplicate entries (listing the same cell
  # twice), we remove these and report.
  # Please check how your csv file is separated and adjust the sep argument
  # in the function if necessary. In this example we import a semicolon
  # separated file.
  read.flow_csv <- function(pathIN){
    raw <- read.csv(pathIN, sep=",", header=TRUE, stringsAsFactors=FALSE)
    IND <- which(duplicated(raw))
    # Check for duplicates and report if found:
    if(any(duplicated(raw))){
      cat(paste0("=== Duplicate entries removed in [",pathIN,"]: ",
                 length(IND)," ===\n"))
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
    # Compute required metadata - column names with description -
    # ranges, min, and max settings
    meta <- data.frame(name=dimnames(dat)[[2]],
                       desc=paste(dimnames(dat)[[2]]),
                       range =(apply(apply(dat,2,range),2,diff)),
                       minRange = apply(dat,2,min),
                       maxRange = apply(dat,2,max))
    # Create flowframe
    flowframef <- new("flowFrame",exprs=as.matrix(dat),
                      parameters=AnnotatedDataFrame(meta))
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

  ## Each parameter of interest needs to be arcsinh transformed with an
  # individual cofactor. The cofactor can be deduced from the size of the
  # linear region around zero on a biexponential scale, as plotted in a
  # histogram (in conventional gating software).
  # Choose manual transformation or automated transformation (we prefer manual)
  # Define parameters and cofactors for transformations:
  dir.create("figures", showWarnings = FALSE)
  gc()
  jpeg(file = "figures/automatedcofactors.jpeg")
  automatedcofactors <- estParamFlowVS(dfs_fs, columnNames)
  dev.off()
  try(capture.output(automatedcofactors,
                     file = "dataPPOutput/automatedcofactors.txt"))

  gc()

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels=columnNames,
                               cofactor=automatedcofactors)
  gc()
  rm(dfs_fs)
  rm(automatedcofactors)


  ##############################
  ######## NORMALIZATION #######
  ##############################

  ## To correct for technical inter-sample variation we apply normalization by
  # fdaNorm (which automatically detects the number of peaks)
  # We continue with the manual transformed dataset
  # Select the markers which require normalization (based on the
  # densityplots you generated above). Be aware that you don't remove
  # biological variation!
  gc()
  dfs_fs_t_auto_normfda <- warpSet(dfs_fs_t_auto, stains=columnNames)
  gc()

  ##############################
  ####### EXPORT TO FCS ########
  ##############################

  ## The flowset (dfs_fs_t_auto_normfda) can be exported to individual
  # fcs files

  #Save flowframes wihtin flowset as fcs files using the flowCore package
  write.flowSet(dfs_fs_t_auto_normfda, outdir='dataPPOutput',
                filename = paste0(gsub(".csv", ".fcs",
                                       sampleNames(dfs_fs_t_auto_normfda))))
  gc()

  ##############################
  ########### Plots ############
  ##############################

  # Pre-Normalized Plots
  flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

  figureDirectory <- paste0(getwd(),"/figures/")

  for (columnName in columnNames) {
    gc()
    columnNameFormula <- as.formula(paste(" ~ ", columnName))
    gc()

    gc()
    jpeg(file = paste0(figureDirectory,"densityPlot",
                       str_replace_all(columnName,"\\.",""), ".jpeg"))
    plot <- densityplot(columnNameFormula, dfs_fs_t_auto, main="auto")
    try(print(plot))
    dev.off()
    gc()

    gc()
    jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFSCA.jpeg",
                       str_replace_all(columnName,"\\.",""), ".jpeg"))
    plot <- densityplot(columnNameFormula, dfs_fs_t_auto_normfda, main="auto")
    try(print(plot))
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

convertToDataFrame <- function(directoryName, columnNames, test) {
  workingDirectory <- getwd()

  clinicalData <- read.csv('data/metadata/metadata.csv')

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  ## Optional: when loading clustered fcs files from cytosplore,
  #provide the directory of the text file 'CSPLR_ST.txt'. Cytosplore exports
  # this file upon running the HSNE. This file contains the decoding of the
  # sample numbers.
  pathST <- "X:/Users/guypw/OneDrive/Documents.txt"

  ## Defining a function to read multiple fcs files from a directory 'dir'
  # into a single data.frame:
  # NB: The column in the output named 'fileName' tracks the original file
  # where each cell came from.
  # Optionally perform remapping of column 'CSPLR_ST' holding cytosplore
  # sample numbers to actual names:
  read.flowdat <- function(dir,path_CSPLR_ST=""){
    # Read:
    filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
    flowset <- read.flowSet(files=filepaths, transformation=FALSE,
                            truncate_max_range = FALSE)
    # Transform to data frame:
    x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
    # Map column 'Original' to filename (in this case holding clusters of
    # HSNE):
    filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs",
                                             full.names=FALSE))
    names(filenames) <- sort(unique(x$Original))
    x$fileName <- filenames[as.character(x$Original)]
    # Remove column 'Original':
    x <- x[,-which(colnames(x)=="Original")]
    # Optionally remap Cytosplore sample tags to original filename:
    if(file.exists(path_CSPLR_ST)){
      # Read:
      sampID <- gsub(
        ".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),
                                  function(x) x[1])))
      names(sampID) <- sapply(
        strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
      x$sampleID <- sampID[as.character(x$CSPLR_ST)]
    }
    return(x)
  }

  ## Read fcs files
  # In our example we will read the data which were clustered in Cytosplore
  # (each fcs file is 1 cluster)
  df <- read.flowdat(dir=dirFCS[1],path_CSPLR_ST = pathST)
  gc()

  if (test) {
    df <- df[seq_len(nrow(df)/50),]
  }

  write.csv(df, 'dataPPOutput/rawDf.csv', row.names = FALSE)
  gc()

  updatedColumnNames <- append(columnNames,"fileName")

  df <- df[,updatedColumnNames]
  gc()
  write.csv(df, 'dataPPOutput/columnsOfInterestDf.csv', row.names = FALSE)
  gc()

  df <- tryCatch({merge(df, clinicalData, by.x = "fileName", by.y = "ï..sample_id")},
                 error=function(x) {
                   merge(df, clinicalData, by.x = "fileName", by.y = "sample_id")
                 })
  df["caseControl"][df["caseControl"] == "Case"] <- 1
  df["caseControl"][df["caseControl"] == "Control"] <- 0

  df["fastSlow"][df["fastSlow"] == "Fast"] <- 1
  df["fastSlow"][df["fastSlow"] == "Slow"] <- 0
  df["fastSlow"][df["fastSlow"] == "N/A"] <- -1
  gc()
  write.csv(df, 'dataPPOutput/columnsOfInterestPlusClinicalDataDf.csv', row.names = FALSE)
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
  caseModelFormula <- as.formula(
    paste("caseControl", paste(columnNames, collapse=" + "), sep=" ~ "))

  progressionModelFormula <- as.formula(
    paste("fastSlow",paste(columnNames, collapse=" + "), sep=" ~ "))

  caseModel <- lm(caseModelFormula, data = df)

  print(summary(caseModel))
  try(capture.output(summary(caseModel),
                     file = "dataPPOutput/summaryCaseModel.txt"))


  print(summary(caseModel)$coefficient)
  try(capture.output(summary(caseModel)$coefficient,
                     file = "dataPPOutput/coefficientCaseModel.txt"))


  print(confint(caseModel))
  try(capture.output(confint(caseModel),
                     file = "dataPPOutput/confintCaseModel.txt"))

  fastModel <- lm(progressionModelFormula, data = df[df["fastSlow"] != -1,])


  print(summary(fastModel))
  try(capture.output(summary(fastModel),
                     file = "dataPPOutput/summaryFastModel.txt"))

  print(summary(fastModel)$coefficient)
  try(capture.output(summary(fastModel)$coefficient,
                     file = "dataPPOutput/coefficientFastModel.txt"))

  print(confint(fastModel))
  try(capture.output(confint(fastModel),
                     file = "dataPPOutput/confintFastModel.txt"))

  rm(fastModel)
  rm(caseModel)

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

flowsomClustering <- function(directoryName, columnNames, numberOfClusters,
                              test) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  ## Optional: when loading clustered fcs files from cytosplore,
  #provide the directory of the text file 'CSPLR_ST.txt'. Cytosplore exports
  # this file upon running the HSNE. This file contains the decoding of the
  # sample numbers.
  pathST <- "X:/Users/guypw/OneDrive/Documents.txt"

  ## Defining a function to read multiple fcs files from a directory 'dir'
  # into a single data.frame:
  # NB: The column in the output named 'fileName' tracks the original file
  # where each cell came from.
  # Optionally perform remapping of column 'CSPLR_ST' holding cytosplore
  # sample numbers to actual names:
  read.flowdat <- function(dir,path_CSPLR_ST=""){
    # Read:
    filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
    flowset <- read.flowSet(files=filepaths[1:2], transformation=FALSE,
                            truncate_max_range = FALSE)
    # Transform to data frame:
    x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
    # Map column 'Original' to filename (in this case holding clusters of
    # HSNE):
    filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs",
                                             full.names=FALSE))[1:2]
    names(filenames) <- sort(unique(x$Original))
    x$fileName <- filenames[as.character(x$Original)]
    # Remove column 'Original':
    x <- x[,-which(colnames(x)=="Original")]
    # Optionally remap Cytosplore sample tags to original filename:
    if(file.exists(path_CSPLR_ST)){
      # Read:
      sampID <- gsub(
        ".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),
                                  function(x) x[1])))
      names(sampID) <- sapply(
        strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
      x$sampleID <- sampID[as.character(x$CSPLR_ST)]
    }
    return(x)
  }

  ## Read fcs files
  # In our example we will read the data which were clustered in Cytosplore
  # (each fcs file is 1 cluster)
  fcsDf <- read.flowdat(dir=dirFCS[1],path_CSPLR_ST = pathST)
  gc()

  dir.create("clusteringOutput", showWarnings = FALSE)

  df <- read.csv('dataPPOutput/columnsOfInterestDf.csv')

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  columnIndexes <- c()

  for (columnName in columnNames) {
    index <- grep(columnName, colnames(fcsDf))
    columnIndexes <- append(columnIndexes, index)
  }

  #run flowsom
  flowsom <- FlowSOM(input = dirFCS,
                     transform = FALSE,
                     scale = FALSE,
                     colsToUse = columnIndexes, #provide the columns for the
                     # clustering
                     nClus = numberOfClusters,
                     seed = 100)

  # Get metaclustering per cell
  clusters_flowsom <- as.factor(flowsom$map$mapping[,1])
  levels(clusters_flowsom) <- flowsom$metaclustering

  if (test) {
    clusters_flowsom <- clusters_flowsom[seq_len(nrow(df))]
  }

  #add flowsom clusters to dataframe
  df <- cbind(df, clusters_flowsom)

  write.csv(df, 'clusteringOutput/flowSomDf.csv', row.names = FALSE)
  try(saveRDS(flowsom, file = "clusteringOutput/flowSom.rds"))
  FlowSOMmary(flowsom, plotFile = "clusteringOutput/FlowSOMmary.pdf")
  rm(flowsom)
  rm(clusters_flowsom)
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

phenographClustering <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("clusteringOutput", showWarnings = FALSE)

  df <- read.csv('clusteringOutput/flowSomDf.csv')

  gc()

  phenograph <- Rphenograph(df[,columnNames], k=knn)
  clusters_phenograph <- as.factor(phenograph$membership)

  #add phenograph clusters to expression data frame
  df <- cbind(df, clusters_phenograph)

  write.csv(df, 'clusteringOutput/phenographDf.csv', row.names = FALSE)
  try(saveRDS(phenograph, file = "clusteringOutput/phenograph.rds"))
  rm(phenograph)
  rm(clusters_phenograph)
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

fastPGClustering <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("clusteringOutput", showWarnings = FALSE)

  df <- read.csv('clusteringOutput/phenographDf.csv')

  gc()

  fastPGResults <- FastPG::fastCluster(as.matrix(df[,columnNames]), knn, 4)
  clusters_fast_pg <- as.factor(fastPGResults$communities)


  #add clusters to expression data frame
  df <- cbind(df, clusters_fast_pg)
  colnames(df)

  write.csv(df, 'clusteringOutput/fastPGDf.csv', row.names = FALSE)
  try(saveRDS(fastPGResults, file = "clusteringOutput/fastPGResults.rds"))
  rm(fastPGResults)
  rm(clusters_fast_pg)
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

umapDimReduction <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <- read.csv('clusteringOutput/fastPGDf.csv')

  umap <- umap(df[,columnNames], n_neighbors = knn, min_dist=0.001,
               verbose=TRUE)
  umap<- as.data.frame(umap)
  colnames(umap) <- c('umap_1', 'umap_2')
  df <- cbind(df,umap)

  write.csv(df, 'clusteringOutput/umapDf.csv', row.names = FALSE)
  rm(umap)
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

visuliseUmap <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  figureDirectory <- paste0(getwd(),"/figures/")

  df <- read.csv('clusteringOutput/umapDf.csv')

  viz.umap <- function(dat,param.name,limits=NULL){
    ColVal <- dat[,param.name]
    if(is.null(limits)){
      Lim <- quantile(ColVal,probs=seq(0,1,0.01))[c(2,100)]
      p <- ggplot(dat, aes(x = umap_1, y =umap_2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=Lim, oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
    } else {
      p <- ggplot(dat, aes(x = umap_1, y = umap_2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=c(limits[1],limits[2]), oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
    }
    p
  }

  for (columnName in columnNames) {
    gc()
    jpeg(file = paste0(figureDirectory,"umap",
                       str_replace_all(columnName,"\\.",""),".jpeg"))
    plot <- viz.umap(dat=df,param.name=columnName)
    try(print(plot))
    dev.off()
    gc()
  }

  #visualize and label clusters on umap
  gc()
  label_flowsom_umap <- df%>%group_by(clusters_flowsom)%>%
    select(umap_1, umap_2)%>%summarize_all(mean)
  label_pheno_umap <- df%>%group_by(clusters_phenograph)%>%
    select(umap_1, umap_2)%>%summarize_all(mean)
  label_fastpg_umap <- df%>%group_by(clusters_fast_pg)%>%
    select(umap_1, umap_2)%>%summarize_all(mean)

  gc()
  jpeg(file = paste0(figureDirectory,"umapFlowsom.jpeg"))
  plot <- ggplot(df, aes(x=umap_1, y=umap_2,color=as.factor(clusters_flowsom)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_flowsom), data=label_flowsom_umap)#+guides(colour=FALSE)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory,"umapPhenograph.jpeg"))
  plot <- ggplot(df, aes(x=umap_1, y=umap_2,color=as.factor(clusters_phenograph)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_phenograph), data=label_pheno_umap)#+guides(colour=FALSE)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory,"umapFastPG.jpeg"))
  plot <- ggplot(df, aes(x=umap_1, y=umap_2,color=as.factor(clusters_fast_pg)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_fast_pg), data=label_fastpg_umap)#+guides(colour=FALSE)
  try(print(plot))
  dev.off()


  rm(label_flowsom_umap)
  rm(label_pheno_umap)
  gc()


  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

diffusionMapDimReduction <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <- read.csv('clusteringOutput/umapDf.csv')

  try(rm(list = c("caseModel","caseModelCoefficient","clinicalData",
                  "fastModel","fastPGResults", "plot")))
  gc()

  # reduce the K, if computational load is too high [it takes approximately 2 hours for the example dataset of 275856 cells]
  dm <- DiffusionMap(df, vars = colnames(df[,columnNames]), k=knn,
                     suppress_dpt = TRUE, verbose=TRUE)

  # add the diffusion components to the expression data frame (either all by dm@eigenvectors, or a selection by dm$DC1, dm$DC2, etc.)
  df <- cbind(df, DC1=dm$DC1, DC2=dm$DC2, DC3=dm$DC3)

  write.csv(df, 'clusteringOutput/diffusionMapDf.csv', row.names = FALSE)
  rm(dm)
  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}

visuliseDiffusionMap <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  figureDirectory <- paste0(getwd(),"/figures/")

  df <- read.csv('clusteringOutput/diffusionMapDf.csv')

  gc()

  viz.dm <- function(dat,dr, param.name,limits=NULL){
    ColVal <- dat[,param.name]
    if(is.null(limits)){
      Lim <- quantile(ColVal,probs=seq(0,1,0.01))[c(2,100)]
      p <- ggplot(dat, aes(x = DC1, y =DC2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=Lim, oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
    } else {
      p <- ggplot(dat, aes(x = DC1, y = DC2)) +geom_point(aes(color = ColVal), size=0.1)+theme_classic()+scale_color_distiller(name=param.name, palette = "RdYlBu", limits=c(limits[1],limits[2]), oob=squish)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ggtitle(param.name)
    }
    p
  }

  for (columnName in columnNames) {
    gc()
    jpeg(file = paste0(figureDirectory,"diffusionMap",
                       str_replace_all(columnName,"\\.",""),".jpeg"))
    plot <- viz.dm(dat=df,param.name=columnName)
    try(print(plot))
    dev.off()
    gc()
  }

  #visualize and label clusters on diffusion map
  #first we need to determine the positions of the cluster labels, based on the DC coordinates
  label_flowsom_dm <- df%>%group_by(clusters_flowsom)%>%
    select(DC1, DC2)%>%summarize_all(mean)
  label_pheno_dm <- df%>%group_by(clusters_phenograph)%>%
    select(DC1, DC2)%>%summarize_all(mean)
  label_fastpg_dm <- df%>%group_by(clusters_fast_pg)%>%
    select(DC1, DC2)%>%summarize_all(mean)

  figureDirectory <- paste0(getwd(),"/figures/")

  gc()
  jpeg(file = paste0(figureDirectory,"diffusionMapFlowsom.jpeg"))
  plot <- ggplot(df, aes(x=DC1, y=DC2, color=as.factor(clusters_flowsom)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_flowsom), data=label_flowsom_dm)#+guides(colour=FALSE)
  print(plot)
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory,"diffusionMapPhenograph.jpeg"))
  plot <- ggplot(df, aes(x=DC1, y=DC2, color=as.factor(clusters_phenograph)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_phenograph), data=label_pheno_dm)#+guides(colour=FALSE)
  print(plot)
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory,"diffusionMapFastPG.jpeg"))
  plot <- ggplot(df, aes(x=DC1, y=DC2, color=as.factor(clusters_fast_pg)))+geom_point(size=0.1)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+geom_label_repel(aes(label=clusters_fast_pg), data=label_fastpg_dm)#+guides(colour=FALSE)
  print(plot)
  dev.off()
  gc()

  rm(list = c("label_flowsom_dm","label_pheno_dm","label_fastpg_dm","plot"))

  gc()

  tryCatch({
    setwd(workingDirectory)},
    error=function(cond) {
      setwd("..")
      setwd("..")
    })
}
