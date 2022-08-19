installlibraries <- function() {
  ## INSTALL REQUIRED PACKAGES ##
  ###############################
  install.packages(
    c(
      "tidyverse",
      "stringr",
      "stringi",
      "knitr",
      "roxygen2",
      "BiocManager",
      "dplyr",
      "R.utils",
      "reshape2",
      "ggplot2",
      "uwot",
      "ggrepel",
      "dplyr",
      "ggplot2",
      "scales",
      "reshape2",
      "RColorBrewer",
      "devtools"
    )
  )
  install.packages(
    c(
      "stringr",
      "locfit",
      "hdrcde",
      "rainbow",
      "fds",
      "fda",
      "flowStats",
      "openCyto",
      "CytoML"
    )
  )

  library(withr)
  library(devtools)
  with_libpaths(new = "R/x86_64-pc-linux-gnu-library/4.2/", install_github("igraph/rigraph@master"))
  with_libpaths(new = "R/x86_64-pc-linux-gnu-library/4.2/", install_github("gastonstat/colortools"))

  BiocManager::install(
    c(
      "flowStats",
      "openCyto",
      "CytoML",
      "Biobase",
      "flowCore",
      "flowVS",
      "flowStats",
      "FlowSOM",
      "slingshot",
      "flowCore",
      "SingleCellExperiment",
      "diffcyt",
      "ComplexHeatmap"
    )
  )
  remotes::install_github("igraph/rigraph@master-old")

  library(devtools)
  devtools::install_github("gastonstat/colortools")
  devtools::install_github("JinmiaoChenLab/cytofkit2")
  devtools::install_github('flying-sheep/knn.covertree')
  devtools::install_github('theislab/destiny')
  devtools::install_github('sararselitsky/FastPG')

  install.packages(c("factoextra", "NbClust", "apcluster"))
  install.packages(c(
    "foreach",
    "doParallel",
    "ranger",
    "palmerpenguins",
    "kableExtra"
  ))
  install.packages("doSNOW")
  install.packages("doParallel")
  install.packages("doMPI")
  install.packages("pheatmap")

  gc()
}

loadlibraries <- function() {
  try(library(ComplexHeatmap))
  try(library(igraph, lib.loc = "R/x86_64-pc-linux-gnu-library/4.2/"))
  try(library(colortools, lib.loc = "R/x86_64-pc-linux-gnu-library/4.2/"))
  try(library(tibble))
  try(library(pheatmap))
  try(library(foreach))
  try(library(doParallel))
  try(library(ranger))
  try(library(palmerpenguins))
  try(library(kableExtra))
  try(library(parallel))
  try(library(R.utils))
  try(library(stringr))
  try(library(flowCore))
  try(library(Biobase))
  try(library(dplyr))
  try(library(flowVS))
  try(library(flowStats))
  try(library(R.utils))
  try(library(flowCore))
  try(library(FlowSOM))
  try(library(SingleCellExperiment))
  try(library(dplyr))
  try(library(ggplot2))
  try(library(scales))
  try(library(reshape2))
  try(library(RColorBrewer))
  try(library(destiny))
  try(library(uwot))
  try(library(slingshot))
  try(library(cytofkit2))
  try(library(ggrepel))
  try(library(tidyverse))
  try(library(diffcyt))
  try(library(corrplot))
  try(library(cluster))
  try(library(factoextra))
  try(library(NbClust))
  try(require(vegan))
  try(library(apcluster))
  try(library(mclust))
}

ungzipFiles <- function() {
  workingDirectory <- getwd()

  dataDirectorys <- c("/data/bCells",
                      "/data/monocytes",
                      "/data/senescence",
                      "/data/tCells")

  for (directory in dataDirectorys) {
    try(setwd(paste0(workingDirectory, directory)))

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

  dataDirectorys <- c("/data/bCells",
                      "/data/monocytes",
                      "/data/senescence",
                      "/data/tCells")

  for (directory in dataDirectorys) {
    try(setwd(paste0(workingDirectory, directory)))

    filenames <- list.files(pattern = ".csv")

    for (filename in filenames) {
      gzip(filename)
    }
  }

  setwd(workingDirectory)

  gc()
}

preprocessing <- function(directoryName,
                          columnNames,
                          test,
                          gate = FALSE,
                          gateTogether = FALSE,
                          gateColumns = NULL,
                          automatedcofactors = NULL) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  # Create an 'output' folder
  gc()
  dir.create("dataPPOutput", showWarnings = FALSE)
  gc()

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  if (test) {
    filenames <- filenames[1:4]
  }

  ## Defining a function to read a flow cytrometry file in csv format:
  # Each row is a cell, each column is a parameter. In our experience, the
  # flow cytometers sometimes output duplicate entries (listing the same cell
  # twice), we remove these and report.
  # Please check how your csv file is separated and adjust the sep argument
  # in the function if necessary. In this example we import a semicolon
  # separated file.
  read.flow_csv <- function(pathIN) {
    raw <-
      read.csv(
        pathIN,
        sep = ",",
        header = TRUE,
        stringsAsFactors = FALSE
      )
    IND <- which(duplicated(raw))
    # Check for duplicates and report if found:
    if (any(duplicated(raw))) {
      cat(paste0(
        "=== Duplicate entries removed in [",
        pathIN,
        "]: ",
        length(IND),
        " ===\n"
      ))
      print(head(raw[IND, ]))
      cat("----\n")
    }
    return(unique(raw))
  }

  # Read all:
  dfs <- sapply(filenames, read.flow_csv, simplify = FALSE)
  gc()

  ##############################
  #REWRITE TO FLOWFRAME/FLOWSET#
  ##############################

  ## Defining a function to rewrite a csv into a flowframe:
  csv_2_ff <- function(dat) {
    # Compute required metadata - column names with description -
    # ranges, min, and max settings
    meta <- data.frame(
      name = dimnames(dat)[[2]],
      desc = paste(dimnames(dat)[[2]]),
      range = (apply(apply(dat, 2, range), 2, diff)),
      minRange = apply(dat, 2, min),
      maxRange = apply(dat, 2, max)
    )
    # Create flowframe
    flowframef <- new("flowFrame",
                      exprs = as.matrix(dat),
                      parameters = AnnotatedDataFrame(meta))
    return(flowframef)
  }

  # rewrite to flowframe
  dfs_ff = sapply(dfs, function(x)
    csv_2_ff(x), simplify = FALSE)
  gc()
  rm(dfs)

  # rewrite to flowset
  dfs_fs <- as(dfs_ff, "flowSet")
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
  if (is.null(automatedcofactors)) {
    dir.create("figures", showWarnings = FALSE)
    gc()
    jpeg(file = "figures/automatedcofactors.jpeg")
    automatedcofactors <- estParamFlowVS(dfs_fs, columnNames)
    dev.off()
    try(capture.output(automatedcofactors,
                       file = "dataPPOutput/automatedcofactors.txt"))

    gc()
  }

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels = columnNames,
                               cofactor = automatedcofactors)
  gc()
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
  dfs_fs_t_auto_normfda <-
    warpSet(dfs_fs_t_auto, stains = columnNames)
  gc()

  # Gate
  if (gate) {
    if (gateTogether) {
      dfs_fs_t_auto_normfda_gated <-
        fsApply(dfs_fs_t_auto_normfda,
                gateTwoMarkersCombinedFcs,
                gateColumns)
    } else {
      dfs_fs_t_auto_normfda_gated <-
        fsApply(dfs_fs_t_auto_normfda, gateMarkersFcs, gateColumns)
    }
  } else {
    dfs_fs_t_auto_normfda_gated <- dfs_fs_t_auto_normfda
  }

  ##############################
  ####### EXPORT TO FCS ########
  ##############################

  ## The flowset (dfs_fs_t_auto_normfda) can be exported to individual
  # fcs files

  #Save flowframes wihtin flowset as fcs files using the flowCore package
  write.flowSet(dfs_fs_t_auto_normfda_gated,
                outdir = 'dataPPOutput',
                filename = paste0(gsub(
                  ".csv", ".fcs",
                  sampleNames(dfs_fs_t_auto_normfda)
                )))

  ##############################
  ########### Plots ############
  ##############################

  # Pre-Normalized Plots
  flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

  figureDirectory <- paste0(getwd(), "/figures/")

  # columnName <- columnNames[1]

  for (columnName in columnNames) {
    gc()
    columnNameFormula <- as.formula(paste(" ~ ", columnName))
    gc()

    gc()
    jpeg(file = paste0(
      figureDirectory,
      "rawDensityPlot",
      str_replace_all(columnName, "\\.", ""),
      ".jpeg"
    ))
    plot <-
      densityplot(columnNameFormula, dfs_fs, main = "auto")
    try(print(plot))
    dev.off()
    gc()

    gc()
    jpeg(file = paste0(
      figureDirectory,
      "transformedDensityPlot",
      str_replace_all(columnName, "\\.", ""),
      ".jpeg"
    ))
    plot <-
      densityplot(columnNameFormula, dfs_fs_t_auto, main = "auto")
    try(print(plot))
    dev.off()
    gc()

    gc()
    jpeg(
      file = paste0(
        figureDirectory,
        "normalisedTransformedDensityPlot",
        str_replace_all(columnName, "\\.", ""),
        ".jpeg"
      )
    )
    plot <-
      densityplot(columnNameFormula, dfs_fs_t_auto_normfda, main = "auto")
    try(print(plot))
    dev.off()
    gc()

    gc()
    jpeg(
      file = paste0(
        figureDirectory,
        "gatedNormalisedTransformedDensityPlot",
        str_replace_all(columnName, "\\.", ""),
        ".jpeg"
      )
    )
    plot <-
      densityplot(columnNameFormula, dfs_fs_t_auto_normfda_gated, main = "auto")
    try(print(plot))
    dev.off()
    gc()
  }

  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

convertToDataFrame <- function(directoryName, columnNames, test) {
  workingDirectory <- getwd()

  clinicalData <- read.csv('data/metadata/metadata.csv')
  head(clinicalData)

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
  read.flowdat <- function(dir, path_CSPLR_ST = "") {
    # Read:
    filepaths <-
      list.files(path = dir,
                 pattern = ".fcs",
                 full.names = TRUE)
    flowset <- read.flowSet(
      files = filepaths,
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    # Transform to data frame:
    x <-
      as.data.frame(exprs(as(flowset, 'flowFrame')), stringsAsFactors = FALSE)
    # Map column 'Original' to filename (in this case holding clusters of
    # HSNE):
    filenames <-
      gsub("[.fcs]",
           "",
           list.files(
             path = dir,
             pattern = ".fcs",
             full.names = FALSE
           ))
    names(filenames) <- sort(unique(x$Original))
    x$fileName <- filenames[as.character(x$Original)]
    # Remove column 'Original':
    x <- x[, -which(colnames(x) == "Original")]
    # Optionally remap Cytosplore sample tags to original filename:
    if (file.exists(path_CSPLR_ST)) {
      # Read:
      sampID <-
        gsub(".fcs", "", basename(sapply(strsplit(readLines(path_CSPLR_ST), ": "),
                                         function(x)
                                           x[1])))
      names(sampID) <-
        sapply(strsplit(readLines(path_CSPLR_ST), ": "), function(x)
          x[2])
      x$sampleID <- sampID[as.character(x$CSPLR_ST)]
    }
    return(x)
  }

  ## Read fcs files
  # In our example we will read the data which were clustered in Cytosplore
  # (each fcs file is 1 cluster)
  df <- read.flowdat(dir = dirFCS[1], path_CSPLR_ST = pathST)
  gc()

  if (test) {
    df <- df[seq_len(nrow(df) / 50), ]
  }

  write.csv(df, 'dataPPOutput/rawDf.csv', row.names = FALSE)
  gc()

  updatedColumnNames <- append(columnNames, "fileName")

  df <- df[, updatedColumnNames]
  gc()
  write.csv(df, 'dataPPOutput/columnsOfInterestDf.csv', row.names = FALSE)
  gc()

  df <- tryCatch({
    merge(df, clinicalData, by.x = "fileName",
          by.y = "Ã¯..patient_id")
  }, error = function(x) {
    merge(df, clinicalData,
          by.x = "fileName", by.y = "patient_id")
  })
  df["caseControl"][df["caseControl"] == "Case"] <- 1
  df["caseControl"][df["caseControl"] == "Control"] <- 0

  df["fastSlow"][df["fastSlow"] == "Fast"] <- 1
  df["fastSlow"][df["fastSlow"] == "Slow"] <- 0
  df["fastSlow"][df["fastSlow"] == "N/A"] <- -1
  gc()
  write.csv(df,
            'dataPPOutput/columnsOfInterestPlusClinicalDataDf.csv',
            row.names = FALSE)
  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}


multipleRegressionTesting <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <-
    read.csv('dataPPOutput/columnsOfInterestPlusClinicalDataDf.csv')

  # This returns the formula:
  caseModelFormula <-
    as.formula(paste("caseControl", paste(columnNames, collapse =
                                            " + "), sep = " ~ "))

  progressionModelFormula <-
    as.formula(paste("fastSlow", paste(columnNames, collapse =
                                         " + "), sep = " ~ "))

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

  fastModel <-
    lm(progressionModelFormula, data = df[df["fastSlow"] != -1, ])


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
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

flowsomClustering <-
  function(directoryName,
           columnNames,
           numberOfClusters,
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
    read.flowdat <- function(dir, path_CSPLR_ST = "") {
      # Read:
      filepaths <-
        list.files(path = dir,
                   pattern = ".fcs",
                   full.names = TRUE)
      flowset <-
        read.flowSet(
          files = filepaths[1:2],
          transformation = FALSE,
          truncate_max_range = FALSE
        )
      # Transform to data frame:
      x <-
        as.data.frame(exprs(as(flowset, 'flowFrame')), stringsAsFactors = FALSE)
      # Map column 'Original' to filename (in this case holding clusters of
      # HSNE):
      filenames <-
        gsub("[.fcs]",
             "",
             list.files(
               path = dir,
               pattern = ".fcs",
               full.names = FALSE
             ))[1:2]
      names(filenames) <- sort(unique(x$Original))
      x$fileName <- filenames[as.character(x$Original)]
      # Remove column 'Original':
      x <- x[, -which(colnames(x) == "Original")]
      # Optionally remap Cytosplore sample tags to original filename:
      if (file.exists(path_CSPLR_ST)) {
        # Read:
        sampID <-
          gsub(".fcs", "", basename(sapply(strsplit(readLines(path_CSPLR_ST), ": "),
                                           function(x)
                                             x[1])))
        names(sampID) <-
          sapply(strsplit(readLines(path_CSPLR_ST), ": "), function(x)
            x[2])
        x$sampleID <- sampID[as.character(x$CSPLR_ST)]
      }
      return(x)
    }

    ## Read fcs files
    # In our example we will read the data which were clustered in Cytosplore
    # (each fcs file is 1 cluster)
    fcsDf <- read.flowdat(dir = dirFCS[1], path_CSPLR_ST = pathST)
    gc()

    dir.create("clusteringOutput", showWarnings = FALSE)

    df <- read.csv('dataPPOutput/columnsOfInterestDf.csv')

    dirFCS <- paste0(getwd(), "/dataPPOutput")

    columnIndexes <- c()

    for (columnName in columnNames) {
      index <- grep(columnName, colnames(fcsDf))
      columnIndexes <- append(columnIndexes, index)
    }

    seed <- 8

    #run flowsom
    flowsom <- FlowSOM(
      input = dirFCS,
      transform = FALSE,
      scale = FALSE,
      colsToUse = columnIndexes,
      #provide the columns for the
      # clustering
      nClus = numberOfClusters,
      seed = seed
    )

    # Get metaclustering per cell
    clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
    meta_clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
    levels(meta_clusters_flowsom) <- flowsom$metaclustering

    if (test) {
      clusters_flowsom <- clusters_flowsom[seq_len(nrow(df))]
      meta_clusters_flowsom <-
        meta_clusters_flowsom[seq_len(nrow(df))]
    }

    #add flowsom clusters to dataframe
    df <- cbind(df, clusters_flowsom)
    df <- cbind(df, meta_clusters_flowsom)

    write.csv(df,
              paste0('clusteringOutput/flowSomDf', numberOfClusters , '.csv'),
              row.names = FALSE)
    try(saveRDS(flowsom,
                file = paste0("clusteringOutput/flowSom", numberOfClusters, ".rds")))
    FlowSOMmary(flowsom,
                plotFile = paste0("clusteringOutput/FlowSOMmary", numberOfClusters, ".pdf"))
    rm(flowsom)
    rm(clusters_flowsom)
    gc()

    tryCatch({
      setwd(workingDirectory)
    },
    error = function(cond) {
      setwd("..")
      setwd("..")
    })
  }

flowsomClusteringFunCluster <-
  function(df, numberOfClusters) {
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
    read.flowdat <- function(dir, path_CSPLR_ST = "") {
      # Read:
      filepaths <-
        list.files(path = dir,
                   pattern = ".fcs",
                   full.names = TRUE)
      flowset <-
        read.flowSet(
          files = filepaths[1:2],
          transformation = FALSE,
          truncate_max_range = FALSE
        )
      # Transform to data frame:
      x <-
        as.data.frame(exprs(as(flowset, 'flowFrame')), stringsAsFactors = FALSE)
      # Map column 'Original' to filename (in this case holding clusters of
      # HSNE):
      filenames <-
        gsub("[.fcs]",
             "",
             list.files(
               path = dir,
               pattern = ".fcs",
               full.names = FALSE
             ))[1:2]
      names(filenames) <- sort(unique(x$Original))
      x$fileName <- filenames[as.character(x$Original)]
      # Remove column 'Original':
      x <- x[, -which(colnames(x) == "Original")]
      # Optionally remap Cytosplore sample tags to original filename:
      if (file.exists(path_CSPLR_ST)) {
        # Read:
        sampID <-
          gsub(".fcs", "", basename(sapply(strsplit(readLines(path_CSPLR_ST), ": "),
                                           function(x)
                                             x[1])))
        names(sampID) <-
          sapply(strsplit(readLines(path_CSPLR_ST), ": "), function(x)
            x[2])
        x$sampleID <- sampID[as.character(x$CSPLR_ST)]
      }
      return(x)
    }

    ## Read fcs files
    # In our example we will read the data which were clustered in Cytosplore
    # (each fcs file is 1 cluster)
    fcsDf <- read.flowdat(dir = dirFCS[1], path_CSPLR_ST = pathST)
    gc()

    dir.create("clusteringOutput", showWarnings = FALSE)

    columnIndexes <- c()

    for (columnName in columnNames) {
      index <- grep(columnName, colnames(fcsDf))
      columnIndexes <- append(columnIndexes, index)
    }

    #run flowsom
    flowsom <- FlowSOM(
      input = dirFCS,
      transform = FALSE,
      scale = FALSE,
      colsToUse = columnIndexes,
      #provide the columns for the
      # clustering
      nClus = numberOfClusters,
      seed = 100
    )

    # Get metaclustering per cell
    cluster <- as.factor(flowsom$map$mapping[, 1])
    levels(cluster) <- flowsom$metaclustering

    results <- list()

    results$cluster <- as.double(cluster)

    return(results)
  }

phenographClustering <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("clusteringOutput", showWarnings = FALSE)

  df <- read.csv('dataPPOutput/columnsOfInterestDf.csv')

  gc()

  phenograph <- Rphenograph(df[, columnNames], k = knn)
  clusters_phenograph <-
    tryCatch({
      as.factor(phenograph$membership)
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      return(phenograph$membership)
    })

  #add phenograph clusters to expression data frame
  df <- cbind(df, clusters_phenograph)

  write.csv(df, 'clusteringOutput/phenographDf.csv', row.names = FALSE)
  try(saveRDS(phenograph, file = "clusteringOutput/phenograph.rds"))
  rm(phenograph)
  rm(clusters_phenograph)
  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

fastPGClustering <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("clusteringOutput", showWarnings = FALSE)

  df <- read.csv('dataPPOutput/columnsOfInterestDf.csv')

  gc()

  fastPGResults <-
    FastPG::fastCluster(as.matrix(df[, columnNames]), knn, 4)
  clusters_fast_pg <-
    tryCatch({
      as.factor(fastPGResults$communities)
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      return(fastPGResults$communities)
    })


  #add clusters to expression data frame
  df <- cbind(df, clusters_fast_pg)

  write.csv(df, 'clusteringOutput/fastPGDf.csv', row.names = FALSE)
  try(saveRDS(fastPGResults, file = "clusteringOutput/fastPGResults.rds"))
  rm(fastPGResults)
  rm(clusters_fast_pg)
  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

umapDimReduction <- function(directoryName, columnNames, knn) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  flowSomDf <-
    read.csv('clusteringOutput/flowSomDf.csv')
  phenographDf <-
    read.csv('clusteringOutput/phenographDf.csv')
  fastPgDf <-
    read.csv('clusteringOutput/fastPGDf.csv')

  df <- flowSomDf
  df[, "clusters_phenograph"] <-
    phenographDf[, "clusters_phenograph"]
  df[, "clusters_fast_pg"] <-
    fastPgDf[, "clusters_fast_pg"]

  randomNumbers <-
    sample(seq(nrow(df)), 50000, replace = FALSE)

  df <- df[randomNumbers,]

  umap <- umap(df[, columnNames],
               n_neighbors = knn,
               min_dist = 0.001,
               verbose = TRUE)
  umap <- as.data.frame(umap)
  colnames(umap) <-
    c('umap_1', 'umap_2')
  df <- cbind(df, umap)

  write.csv(df, 'clusteringOutput/umapDf.csv', row.names = FALSE)
  rm(umap)
  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

visuliseUmap <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  figureDirectory <- paste0(getwd(), "/figures/")

  df <- read.csv('clusteringOutput/umapDf.csv')

  viz.umap <- function(dat, param.name, limits = NULL) {
    ColVal <- dat[, param.name]
    if (is.null(limits)) {
      Lim <- quantile(ColVal, probs = seq(0, 1, 0.01))[c(2, 100)]
      p <-
        ggplot(dat, aes(x = umap_1, y = umap_2)) + geom_point(aes(color = ColVal), size =
                                                                0.1) + theme_classic() + scale_color_distiller(
                                                                  name = param.name,
                                                                  palette = "RdYlBu",
                                                                  limits = Lim,
                                                                  oob = squish
                                                                ) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank()) + ggtitle(param.name)
    } else {
      p <-
        ggplot(dat, aes(x = umap_1, y = umap_2)) + geom_point(aes(color = ColVal), size =
                                                                0.1) + theme_classic() + scale_color_distiller(
                                                                  name = param.name,
                                                                  palette = "RdYlBu",
                                                                  limits = c(limits[1], limits[2]),
                                                                  oob = squish
                                                                ) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank()) + ggtitle(param.name)
    }
    p
  }

  for (columnName in columnNames) {
    gc()
    jpeg(file = paste0(
      figureDirectory,
      "umap",
      str_replace_all(columnName, "\\.", ""),
      ".jpeg"
    ))
    plot <- viz.umap(dat = df, param.name = columnName)
    try(print(plot))
    dev.off()
    gc()
  }

  metaFlowSomMarkers <-
    read.csv('clusteringOutput/meta_clusters_flowsomMarkers.csv')
  colnames(metaFlowSomMarkers)[length(colnames(metaFlowSomMarkers))] <-
    "meta_flowsom_markers"
  df <-
    merge(
      x = df,
      y = metaFlowSomMarkers[, c("meta_clusters_flowsom", "meta_flowsom_markers")],
      by.x = "meta_clusters_flowsom",
      by.y = "meta_clusters_flowsom",
      all.x = TRUE
    )

  flowSomMarkers <-
    read.csv('clusteringOutput/clusters_flowsomMarkers.csv')
  colnames(flowSomMarkers)[length(colnames(flowSomMarkers))] <-
    "flowsom_markers"
  df <-
    merge(
      x = df,
      y = flowSomMarkers[, c("clusters_flowsom", "flowsom_markers")],
      by.x = "clusters_flowsom",
      by.y = "clusters_flowsom",
      all.x = TRUE
    )

  fastPgMarkers <-
    read.csv('clusteringOutput/clusters_fast_pgMarkers.csv')
  colnames(fastPgMarkers)[length(colnames(fastPgMarkers))] <-
    "fastpg_markers"
  df <-
    merge(
      x = df,
      y = fastPgMarkers[, c("clusters_fast_pg", "fastpg_markers")],
      by.x = "clusters_fast_pg",
      by.y = "clusters_fast_pg",
      all.x = TRUE
    )

  phenographMarkers <-
    read.csv('clusteringOutput/clusters_phenographMarkers.csv')
  colnames(phenographMarkers)[length(colnames(phenographMarkers))] <-
    "phenograph_markers"
  df <-
    merge(
      x = df,
      y = phenographMarkers[, c("clusters_phenograph", "phenograph_markers")],
      by.x = "clusters_phenograph",
      by.y = "clusters_phenograph",
      all.x = TRUE
    )


  metaFlowSomcellPopulations <-
    read.csv('clusteringOutput/meta_clusters_flowsomCellPopulations.csv')
  colnames(metaFlowSomcellPopulations)[length(colnames(metaFlowSomcellPopulations))] <-
    "meta_flowsom_cell_population"
  df <-
    merge(
      x = df,
      y = metaFlowSomcellPopulations[, c("meta_clusters_flowsom", "meta_flowsom_cell_population")],
      by.x = "meta_clusters_flowsom",
      by.y = "meta_clusters_flowsom",
      all.x = TRUE
    )

  flowSomcellPopulations <-
    read.csv('clusteringOutput/clusters_flowsomCellPopulations.csv')
  colnames(flowSomcellPopulations)[length(colnames(flowSomcellPopulations))] <-
    "flowsom_cell_population"
  df <-
    merge(
      x = df,
      y = flowSomcellPopulations[, c("clusters_flowsom", "flowsom_cell_population")],
      by.x = "clusters_flowsom",
      by.y = "clusters_flowsom",
      all.x = TRUE
    )

  fastPgcellPopulations <-
    read.csv('clusteringOutput/clusters_fast_pgCellPopulations.csv')
  colnames(fastPgcellPopulations)[length(colnames(fastPgcellPopulations))] <-
    "fastpg_cell_population"
  df <-
    merge(
      x = df,
      y = fastPgcellPopulations[, c("clusters_fast_pg", "fastpg_cell_population")],
      by.x = "clusters_fast_pg",
      by.y = "clusters_fast_pg",
      all.x = TRUE
    )

  phenographcellPopulations <-
    read.csv('clusteringOutput/clusters_phenographCellPopulations.csv')
  colnames(phenographcellPopulations)[length(colnames(phenographcellPopulations))] <-
    "phenograph_cell_population"
  df <-
    merge(
      x = df,
      y = phenographcellPopulations[, c("clusters_phenograph", "phenograph_cell_population")],
      by.x = "clusters_phenograph",
      by.y = "clusters_phenograph",
      all.x = TRUE
    )


  #visualize and label clusters on umap
  gc()
  label_flowsom_umap <- df %>% group_by(flowsom_cell_population) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_flowsom_cell_population) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_pheno_umap <-
    df %>% group_by(phenograph_cell_population) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_fastpg_umap <- df %>% group_by(fastpg_cell_population) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomCellPopulations.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(flowsom_cell_population)
    )) + geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = flowsom_cell_population), data = label_flowsom_umap)

  try(print(plot))
  dev.off()
  gc()
  gc()
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomCellPopulations.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(meta_flowsom_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = meta_flowsom_cell_population), data = label_meta_flowsom_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapPhenographCellPopulations.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(phenograph_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = phenograph_cell_population), data = label_pheno_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapFastPGCellPopulations.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(fastpg_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = fastpg_cell_population), data = label_fastpg_umap)

  try(print(plot))
  dev.off()

  label_flowsom_umap <- df %>% group_by(clusters_flowsom) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_clusters_flowsom) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_pheno_umap <- df %>% group_by(clusters_phenograph) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_fastpg_umap <- df %>% group_by(clusters_fast_pg) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(clusters_flowsom)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_flowsom), data = label_flowsom_umap)

  try(print(plot))
  dev.off()
  gc()
  gc()
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(meta_clusters_flowsom)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = meta_clusters_flowsom), data = label_meta_flowsom_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapPhenographClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(clusters_phenograph)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_phenograph), data = label_pheno_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapFastPGClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(clusters_fast_pg)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_fast_pg), data = label_fastpg_umap)

  try(print(plot))
  dev.off()

  label_flowsom_umap <- df %>% group_by(flowsom_markers) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_flowsom_markers) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_pheno_umap <- df %>% group_by(phenograph_markers) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)
  label_fastpg_umap <- df %>% group_by(fastpg_markers) %>%
    select(umap_1, umap_2) %>% summarize_all(mean)

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomMarkers.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(flowsom_markers)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = flowsom_markers), data = label_flowsom_umap)

  try(print(plot))
  dev.off()
  gc()
  gc()
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomMarkers.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(meta_flowsom_markers)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = meta_flowsom_markers), data = label_meta_flowsom_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapPhenographMarkers.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(phenograph_markers)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = phenograph_markers), data = label_pheno_umap)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "umapFastPGMarkers.jpeg"))
  plot <-
    ggplot(df, aes(
      x = umap_1,
      y = umap_2,
      color = as.factor(fastpg_markers)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = fastpg_markers), data = label_fastpg_umap)

  try(print(plot))
  dev.off()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

diffusionMapDimReduction <-
  function(directoryName, columnNames, knn) {
    workingDirectory <- getwd()

    setwd(paste0("./data/", directoryName))

    flowSomDf <-
      read.csv('clusteringOutput/flowSomDf.csv')
    phenographDf <-
      read.csv('clusteringOutput/phenographDf.csv')
    fastPgDf <-
      read.csv('clusteringOutput/fastPGDf.csv')

    df <- flowSomDf
    df[, "clusters_phenograph"] <-
      phenographDf[, "clusters_phenograph"]
    df[, "clusters_fast_pg"] <-
      fastPgDf[, "clusters_fast_pg"]

    randomNumbers <-
      sample(seq(nrow(df)), 100000, replace = FALSE)

    df <- df[randomNumbers,]

    gc()

    # reduce the K, if computational load is too high [it takes approximately 2 hours for the example dataset of 275856 cells]
    dm <- DiffusionMap(
      df,
      vars = colnames(df[, columnNames]),
      k = knn,
      suppress_dpt = TRUE,
      verbose = TRUE
    )

    # add the diffusion components to the expression data frame (either all by dm@eigenvectors, or a selection by dm$DC1, dm$DC2, etc.)
    df <- cbind(df,
                DC1 = dm$DC1,
                DC2 = dm$DC2,
                DC3 = dm$DC3)

    write.csv(df, 'clusteringOutput/diffusionMapDf.csv', row.names = FALSE)
    rm(dm)
    gc()

    tryCatch({
      setwd(workingDirectory)
    },
    error = function(cond) {
      setwd("..")
      setwd("..")
    })
  }

visuliseDiffusionMap <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  figureDirectory <- paste0(getwd(), "/figures/")

  df <- read.csv('clusteringOutput/diffusionMapDf.csv')

  gc()

  viz.dm <- function(dat, dr, param.name, limits = NULL) {
    ColVal <- dat[, param.name]
    if (is.null(limits)) {
      Lim <- quantile(ColVal, probs = seq(0, 1, 0.01))[c(2, 100)]
      p <-
        ggplot(dat, aes(x = DC1, y = DC2)) + geom_point(aes(color = ColVal), size =
                                                          0.1) + theme_classic() + scale_color_distiller(
                                                            name = param.name,
                                                            palette = "RdYlBu",
                                                            limits = Lim,
                                                            oob = squish
                                                          ) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                 panel.grid.minor = element_blank()) + ggtitle(param.name)
    } else {
      p <-
        ggplot(dat, aes(x = DC1, y = DC2)) + geom_point(aes(color = ColVal), size =
                                                          0.1) + theme_classic() + scale_color_distiller(
                                                            name = param.name,
                                                            palette = "RdYlBu",
                                                            limits = c(limits[1], limits[2]),
                                                            oob = squish
                                                          ) + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                 panel.grid.minor = element_blank()) + ggtitle(param.name)
    }
    p
  }

  for (columnName in columnNames) {
    gc()
    jpeg(file = paste0(
      figureDirectory,
      "diffusionMap",
      str_replace_all(columnName, "\\.", ""),
      ".jpeg"
    ))
    plot <- viz.dm(dat = df, param.name = columnName)
    try(print(plot))
    dev.off()
    gc()
  }

  metaFlowSomcellPopulations <-
    read.csv('clusteringOutput/meta_clusters_flowsomCellPopulations.csv')
  colnames(metaFlowSomcellPopulations)[length(colnames(metaFlowSomcellPopulations))] <-
    "meta_flowsom_cell_population"
  df <-
    merge(
      x = df,
      y = metaFlowSomcellPopulations[, c("meta_clusters_flowsom", "meta_flowsom_cell_population")],
      by.x = "meta_clusters_flowsom",
      by.y = "meta_clusters_flowsom",
      all.x = TRUE
    )

  flowSomcellPopulations <-
    read.csv('clusteringOutput/clusters_flowsomCellPopulations.csv')
  colnames(flowSomcellPopulations)[length(colnames(flowSomcellPopulations))] <-
    "flowsom_cell_population"
  df <-
    merge(
      x = df,
      y = flowSomcellPopulations[, c("clusters_flowsom", "flowsom_cell_population")],
      by.x = "clusters_flowsom",
      by.y = "clusters_flowsom",
      all.x = TRUE
    )

  fastPgcellPopulations <-
    read.csv('clusteringOutput/clusters_fast_pgCellPopulations.csv')
  colnames(fastPgcellPopulations)[length(colnames(fastPgcellPopulations))] <-
    "fastpg_cell_population"
  df <-
    merge(
      x = df,
      y = fastPgcellPopulations[, c("clusters_fast_pg", "fastpg_cell_population")],
      by.x = "clusters_fast_pg",
      by.y = "clusters_fast_pg",
      all.x = TRUE
    )

  phenographcellPopulations <-
    read.csv('clusteringOutput/clusters_phenographCellPopulations.csv')
  colnames(phenographcellPopulations)[length(colnames(phenographcellPopulations))] <-
    "phenograph_cell_population"
  df <-
    merge(
      x = df,
      y = phenographcellPopulations[, c("clusters_phenograph", "phenograph_cell_population")],
      by.x = "clusters_phenograph",
      by.y = "clusters_phenograph",
      all.x = TRUE
    )


  #visualize and label clusters on umap
  gc()
  label_flowsom_dm <- df %>% group_by(flowsom_cell_population) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_meta_flowsom_dm <-
    df %>% group_by(meta_flowsom_cell_population) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_pheno_dm <- df %>% group_by(phenograph_cell_population) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_fastpg_dm <- df %>% group_by(fastpg_cell_population) %>%
    select(DC1, DC2) %>% summarize_all(mean)

  gc()
  jpeg(file = paste0(
    figureDirectory,
    "diffusionMapFlowsomCellPopulations.jpeg"
  ))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(flowsom_cell_population)
    )) + geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = flowsom_cell_population), data = label_flowsom_dm)
  try(print(plot))
  dev.off()
  gc()
  gc()
  jpeg(file = paste0(
    figureDirectory,
    "diffusionMapMetaFlowsomCellPopulations.jpeg"
  ))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(meta_flowsom_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = meta_flowsom_cell_population), data = label_meta_flowsom_dm)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(
    figureDirectory,
    "diffusionMapPhenographCellPopulations.jpeg"
  ))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(phenograph_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = phenograph_cell_population), data = label_pheno_dm)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapFastPGCellPopulations.jpeg"))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(fastpg_cell_population)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = fastpg_cell_population), data = label_fastpg_dm)

  try(print(plot))
  dev.off()

  label_flowsom_dm <- df %>% group_by(clusters_flowsom) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_meta_flowsom_dm <-
    df %>% group_by(meta_clusters_flowsom) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_pheno_dm <- df %>% group_by(clusters_phenograph) %>%
    select(DC1, DC2) %>% summarize_all(mean)
  label_fastpg_dm <- df %>% group_by(clusters_fast_pg) %>%
    select(DC1, DC2) %>% summarize_all(mean)

  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapFlowsomClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(clusters_flowsom)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_flowsom), data = label_flowsom_dm)

  try(print(plot))
  dev.off()
  gc()
  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapMetaFlowsomClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(meta_clusters_flowsom)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = meta_clusters_flowsom), data = label_meta_flowsom_dm)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapPhenographClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(clusters_phenograph)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_phenograph), data = label_pheno_dm)
  try(print(plot))
  dev.off()
  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapFastPGClusters.jpeg"))
  plot <-
    ggplot(df, aes(
      x = DC1,
      y = DC2,
      color = as.factor(clusters_fast_pg)
    )) +
    geom_point(size = 0.1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "None"
    ) +
    geom_label_repel(aes(label = clusters_fast_pg), data = label_fastpg_dm)
  try(print(plot))
  dev.off()

  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

veganOptimalClusters <-
  function(directoryName,
           columnNames,
           minClusters,
           maxClusters) {
    workingDirectory <- getwd()

    setwd(paste0("./data/", directoryName))

    dirFCS <- paste0(getwd(), "/dataPPOutput")

    df <- read.csv("dataPPOutput/columnsOfInterestDf.csv")

    columnIndexes <- c()

    for (columnName in columnNames) {
      index <- grep(columnName, colnames(df))
      columnIndexes <- append(columnIndexes, index)
    }

    minimalDf <- as.matrix(df[#seq_len(nrow(df)/100)
      , columnIndexes])


    figureDirectory <- paste0(getwd(), "/figures/")

    # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
    fit <-
      cascadeKM(scale(minimalDf, center = TRUE, scale = TRUE),
                minClusters,
                maxClusters,
                iter = 100)
    gc()
    jpeg(file = paste0(
      figureDirectory,
      "vegan",
      as.character(minClusters),
      "-",
      as.character(maxClusters),
      ".jpeg"
    ))
    par(mar = c(1, 1, 1, 1))
    p <- plot(fit, sortg = TRUE, grpmts.plot = TRUE)
    print(p)
    dev.off()
    gc()

    calinski.best <- as.numeric(which.max(fit$results[2,]))

    try(capture.output(
      cat(
        "Calinski criterion optimal number of clusters:",
        calinski.best,
        "\n"
      ),
      file = paste0(
        "clusteringOutput/veganOptimumClusters",
        as.character(minClusters),
        "-",
        as.character(maxClusters),
        ".txt"
      )
    ))



    tryCatch({
      setwd(workingDirectory)
    },
    error = function(cond) {
      setwd("..")
      setwd("..")
    })

    gc()
  }

mclustOptimalClusters <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  df <- read.csv("dataPPOutput/columnsOfInterestDf.csv")

  columnIndexes <- c()

  for (columnName in columnNames) {
    index <- grep(columnName, colnames(df))
    columnIndexes <- append(columnIndexes, index)
  }

  minimalDf <- as.matrix(df[, columnIndexes])

  figureDirectory <- paste0(getwd(), "/figures/")

  # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
  d_clust <- Mclust(minimalDf, G = 1:2)
  m.best <- dim(d_clust$z)[2]

  try(capture.output(
    cat("model-based optimal number of clusters:", m.best, "\n"),
    file = paste0(
      "clusteringOutput/mclustOptimumClusters.txt",
      clusterName,
      columnName ,
      ".txt"
    )
  ))

  jpeg(file = paste0(figureDirectory, "mclustBic.jpeg"))
  par(mar = c(1, 1, 1, 1))
  p <- plot(d_clust$BIC)
  print(p)
  dev.off()
  gc()

  # Manual step
  # p <- plot(d_clust)

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })

  gc()
}

apclusterOptimalClusters <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  df <- read.csv("dataPPOutput/columnsOfInterestDf.csv")

  columnIndexes <- c()

  for (columnName in columnNames) {
    index <- grep(columnName, colnames(df))
    columnIndexes <- append(columnIndexes, index)
  }

  minimalDf <- as.matrix(df[, columnIndexes])

  figureDirectory <- paste0(getwd(), "/figures/")

  # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
  d.apclus <- apcluster(negDistMat(r = 2), minimalDf)

  try(capture.output(
    cat(
      "affinity propogation optimal number of clusters:",
      length(d.apclus@clusters),
      "\n"
    ),
    file = paste0(
      "clusteringOutput/apclusterOptimumClusters.txt",
      clusterName,
      columnName ,
      ".txt"
    )
  ))

  # 4
  jpeg(file = paste0(figureDirectory, "apcluster1.jpeg"))
  par(mar = c(1, 1, 1, 1))
  p <- heatmap(d.apclus)
  print(p)
  dev.off()
  gc()

  jpeg(file = paste0(figureDirectory, "apcluster1.jpeg"))
  par(mar = c(1, 1, 1, 1))
  p <- plot(d.apclus, d)
  print(p)
  dev.off()
  gc()

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })

  gc()
}


differentialStatesVolcanoPlot <- function(res_DS_DT,
                                          figureDirectory,
                                          mycolors,
                                          figureTitle) {
  gc()
  jpeg(file = paste0(figureDirectory,
                     figureTitle))
  par(mar = c(1, 1, 1, 1))
  p <- ggplot(data = res_DS_DT,
              aes(
                x = logFC,
                y = minus_log_p_adj,
                col = diff_expressed,
                label = id
              )) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    scale_colour_manual(values = mycolors) +
    geom_text_repel()
  print(p)
  dev.off()
}

differentialAbundanceVolcanoPlot <- function(res_DS_DT,
                                             figureDirectory,
                                             mycolors,
                                             figureTitle) {
  gc()
  jpeg(file = paste0(figureDirectory,
                     figureTitle))
  par(mar = c(1, 1, 1, 1))
  p <- ggplot(data = res_DS_DT,
              aes(
                x = logFC,
                y = minus_log_p_adj,
                col = diff_expressed,
                label = cluster_id
              )) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    scale_colour_manual(values = mycolors) +
    geom_text_repel()
  print(p)
  dev.off()
}

differentialStatesManhattanPlot <- function(res_DS_DT,
                                            figureDirectory,
                                            figureTitle) {
  gc()
  jpeg(file = paste0(figureDirectory,
                     figureTitle))
  par(mar = c(1, 1, 1, 1))
  p <-
    ggplot(res_DS_DT, aes(id, minus_log_p_adj)) +
    geom_point() + geom_hline(yintercept = 0 - log10(0.05)) +
    xlab("Cluster") + ylab("-log10(Adjusted P-Value)") +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
  gc()
}

differentialAbundanceManhattanPlot <- function(res_DA_DT,
                                               figureDirectory,
                                               figureTitle) {
  gc()
  jpeg(file = paste0(figureDirectory,
                     figureTitle))
  par(mar = c(1, 1, 1, 1))
  p <-
    ggplot(res_DA_DT, aes(cluster_id, minus_log_p_adj)) +
    geom_point() + geom_hline(yintercept = 0 - log10(0.05)) +
    xlab("Cluster") + ylab("-log10(Adjusted P-Value)") +
    theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
  gc()
}

differentialStatesSaveResults <- function(res_DS_DT,
                                          res_DS,
                                          res_DA_DT,
                                          res_DA,
                                          datasetTitle) {
  write.csv(
    res_DS_DT,
    paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialStatesStatistics.csv'
    ),
    row.names = FALSE
  )
  try(saveRDS(
    res_DS,
    file = paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialStatesStatistics.rds'
    )
  ))
  write.csv(
    res_DA_DT,
    paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialAbundanceStatistics.csv'
    ),
    row.names = FALSE
  )
  try(saveRDS(
    res_DA,
    file = paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialAbundanceStatistics.rds'
    )
  ))
}

gateMarkers <- function(directoryName,
                        columnNames = NULL,
                        cutoff = NULL) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <- read.csv("dataPPOutput/columnsOfInterestDf.csv")

  if (!is.null(columnNames)) {
    for (columnName in columnNames) {
      df <-
        df[df[, columnName] >= cutoff[match(columnName, columnNames)],]
    }
  }

  write.csv(df, 'dataPPOutput/gatedDf.csv', row.names = FALSE)

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })

  gc()
}

gateTwoMarkersCombined <-
  function(directoryName,
           columnNames = NULL,
           cutoff = NULL) {
    workingDirectory <- getwd()

    setwd(paste0("./data/", directoryName))

    df <- read.csv("dataPPOutput/columnsOfInterestDf.csv")

    if (!is.null(columnNames)) {
      df <-
        df[df[, columnNames[1]] >= cutoff[1] |
             df[, columnNames[2]] >= cutoff[2], ]
    }

    write.csv(df, 'dataPPOutput/gatedDf.csv', row.names = FALSE)

    tryCatch({
      setwd(workingDirectory)
    },
    error = function(cond) {
      setwd("..")
      setwd("..")
    })

    gc()
  }

gateTwoMarkersCombinedFcs <- function(df, gateColumns) {
  columnName1 <- colnames(gateColumns)[1]
  columnName2 <- colnames(gateColumns)[2]

  if (!is.null(gateColumns)) {
    df <-
      df[df[, columnName1] >= gateColumns[, columnName1] |
           df[, columnName2] >= gateColumns[, columnName2], ]
  }

  return(df)

  gc()
}

gateMarkersFcs <- function(df, gateColumns) {
  columnName1 <- colnames(gateColumns)[1]
  columnName2 <- colnames(gateColumns)[2]

  if (!is.null(gateColumns)) {
    df <-
      df[df[, columnName1] >= gateColumns[, columnName1] &
           df[, columnName2] >= gateColumns[, columnName2], ]
  }

  return(df)

  gc()
}

differentialAbundanceAnalysis <- function(directoryName,
                                          columnNames,
                                          clusterName,
                                          samplesContributionToClustersThreshold,
                                          differentialAbundanceThreshold,
                                          calculateSampleContributionsToClusters,
                                          group_id,
                                          visits,
                                          cases,
                                          covariants,
                                          singleCluster,
                                          markersOrCell,
                                          progression = "") {
  concatinatedVisits <- toString(visits)

  # Read experiment data
  experimentInfo <- read.csv("data/metadata/metadata.csv")

  if (markersOrCell != "Clusters") {
    cellPopulationMarkers <-
      read.csv(
        paste0(
          "data/",
          directoryName,
          "/clusteringOutput/",
          clusterName,
          markersOrCell,
          ".csv"
        )
      )
  }

  if (directoryName == "monocytes") {
    experimentInfo[which(experimentInfo[, "sample_id"] == "BAS_057_02", arr.ind =
                           TRUE), "sample_id"] <- "BAS057_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00074-2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00074-02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00075-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00075-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00092-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00092-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_9", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186-09"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-3", arr.ind =
                           TRUE), "sample_id"] <- "BLT00211-03"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00211-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00211-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214-2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00214-02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214-5", arr.ind =
                           TRUE), "sample_id"] <- "BLT00214-05"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00230-2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00230-02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00243-7", arr.ind =
                           TRUE), "sample_id"] <- "BLT00243-07"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00244-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00244-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00244-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00254-5", arr.ind =
                           TRUE), "sample_id"] <- "BLT00254-05"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00254-7", arr.ind =
                           TRUE), "sample_id"] <- "BLT00254-07"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00265-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00265-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00265-8", arr.ind =
                           TRUE), "sample_id"] <- "BLT00265-08"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00271-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00271-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00271-7", arr.ind =
                           TRUE), "sample_id"] <- "BLT00271-07"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00274-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-3", arr.ind =
                           TRUE), "sample_id"] <- "BLT00285-03"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-5", arr.ind =
                           TRUE), "sample_id"] <- "BLT00285-05"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00285-5", arr.ind =
                           TRUE), "sample_id"] <- "BLT00285-05"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT000286-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00286-04"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00286-6", arr.ind =
                           TRUE), "sample_id"] <- "BLT00286-06"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00297_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00297_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00244_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT0244_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "QS_024-2", arr.ind =
                           TRUE), "sample_id"] <- "QS_024-02"
  } else if (directoryName == "senescence") {
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00074_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00074-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00075_06", arr.ind =
                           TRUE), "sample_id"] <- "BLT00075-6"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_9", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186-9"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00261_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00261-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00297_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00297-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00075-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00075-4_R1"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00211-3", arr.ind =
                           TRUE), "sample_id"] <- "BLT00211-3 "
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00297_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00297-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214-5", arr.ind =
                           TRUE), "sample_id"] <- "BLT000214-5"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00057_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00057-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00092_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00092_4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00198_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00198_2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00214_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00214-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00230_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00230-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00242_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00242-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00243_05", arr.ind =
                           TRUE), "sample_id"] <- "BLT00243-5"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00274-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00243_05", arr.ind =
                           TRUE), "sample_id"] <- "BLT00243-5"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00274-2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274-05", arr.ind =
                           TRUE), "sample_id"] <- "BLT00274-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT000286-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00286-5"
  } else if (directoryName == "tCells") {
    experimentInfo[which(experimentInfo[, "sample_id"] == "BAS_057_02", arr.ind =
                           TRUE), "sample_id"] <- "BAS057_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BAS_033_02", arr.ind =
                           TRUE), "sample_id"] <- "BAS033_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BAS_057_02", arr.ind =
                           TRUE), "sample_id"] <- "BAS057_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00074_04", arr.ind =
                           TRUE), "sample_id"] <- "BLT00074-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00075_06", arr.ind =
                           TRUE), "sample_id"] <- "BLT00075-6"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_02", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186_2"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BAS101", arr.ind =
                           TRUE), "sample_id"] <- "BAS00101"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00186_9", arr.ind =
                           TRUE), "sample_id"] <- "BLT00186-9"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00261_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00261_02"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00274-05", arr.ind =
                           TRUE), "sample_id"] <- "BLT00274-5"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT000286-4", arr.ind =
                           TRUE), "sample_id"] <- "BLT00286-4"
    experimentInfo[which(experimentInfo[, "sample_id"] == "BLT00297_2", arr.ind =
                           TRUE), "sample_id"] <- "BLT00297-2"
  }

  experimentInfo <-
    experimentInfo[order(experimentInfo[, "sample_id"]), ]
  experimentInfo <-
    experimentInfo[experimentInfo[, "patient_id"] %in%
                     experimentInfo[experimentInfo[, "visit"] == max(visits), "patient_id"],]
  experimentInfo <-
    experimentInfo[experimentInfo[, "visit"] %in% visits, ]
  experimentInfo <-
    experimentInfo[experimentInfo[, "caseControl"] %in% cases, ]
  if (progression != "") {
    experimentInfo <-
      experimentInfo[experimentInfo[, "fastSlow"] == progression |
                       experimentInfo[, "caseControl"] == "Control", ]
  }
  experimentInfo[, "group_id"] <- NA
  experimentInfo[, "patient_id"] <-
    factor(experimentInfo[, "patient_id"])
  experimentInfo[, "sample_id"] <-
    factor(experimentInfo[, "sample_id"])
  experimentInfo[, "gender"] <-
    factor(experimentInfo[, "gender"])

  experimentInfo <- experimentInfo[!is.na(experimentInfo[, "ethnicity"]),]
  experimentInfo[, "ethnicity"] <-
    factor(experimentInfo[, "ethnicity"])
  experimentInfo[, "ageAtVisitDouble"] <-
    as.double(experimentInfo[, "ageAtVisit"])
  experimentInfo[, "timeFromVisit1InDaysDouble"] <-
    as.double(experimentInfo[, "timeFromVisit1InDays"])
  experimentInfo[, "fastSlow"] <-
    factor(experimentInfo[, "fastSlow"])
  experimentInfo[, "bulbarLimb"] <-
    factor(experimentInfo[, "bulbarLimb"])
  experimentInfo[, "caseControl"] <-
    factor(experimentInfo[, "caseControl"])
  experimentInfo[, "visit"] <-
    as.double(experimentInfo[, "visit"])
  experimentInfo[, "group_id"] <-
    factor(experimentInfo[, group_id])
  experimentInfo[, "patient_id"] <-
    factor(experimentInfo[, "patient_id"])

  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")

  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("differentialTestingOutputs", showWarnings = FALSE)
  figureDirectory <- paste0(getwd(), "/figures/")

  # Read csv
  if (clusterName == "clusters_flowsom") {
    df <- read.csv("clusteringOutput/flowSomDf.csv")
  } else if (clusterName == "meta_clusters_flowsom") {
    df <- read.csv("clusteringOutput/flowSomDf.csv")
  } else if (clusterName == "clusters_phenograph") {
    df <- read.csv("clusteringOutput/phenographDf.csv")
  } else if (clusterName == "clusters_fast_pg") {
    df <- read.csv("clusteringOutput/fastPGDf.csv")
  }

  df <- df[order(df[, "fileName"]), ]

  # Filter df for samples
  df <- df[df[, "fileName"] %in% experimentInfo[, "sample_id"], ]

  if (markersOrCell != "Clusters") {
    df <-
      merge(df,
            cellPopulationMarkers[, c(clusterName, "cell_population")],
            by = clusterName,
            all.x = TRUE)
  }

  # Filter df for samples
  experimentInfo <-
    experimentInfo[experimentInfo[, "sample_id"] %in% df[, "fileName"],]

  # Extract only the relevant columns
  minimalDf <- df[, c(columnNames)]

  # split the dataframe into a dataframe for each file
  listOfDfs <- list()
  for (file in unique(minimalDf[, "fileName"])) {
    minimalDfExtract <- minimalDf[minimalDf[, "fileName"] == file, ]
    minimalDfExtract <-
      minimalDfExtract[,!(names(minimalDfExtract) %in% c("fileName"))]
    listOfDfs <- append(listOfDfs, list(minimalDfExtract))
  }

  # Create marker information
  markerColumnNames <- columnNames[columnNames != "fileName"]
  markerInformation <- data.frame(markerColumnNames)
  markerInformation[, "channel_name"] <- markerColumnNames
  markerInformation[, "marker_name"] <- markerColumnNames
  markerInformation[, "marker_class"] <-
    c(rep("type", length(markerColumnNames) - 2), rep("state", 2))
  markerInformation <-
    markerInformation[, 2:ncol(markerInformation)]

  # Transform the input into the correct format
  d_se <- prepareData(listOfDfs, experimentInfo, markerInformation)
  if (singleCluster) {
    rowData(d_se)[, "cluster_id"] <- 1
    clusterType <- "AllCells"
  } else if (markersOrCell != "Clusters") {
    rowData(d_se)[, "cluster_id"] <- df[, "cell_population"]
    clusterType <- markersOrCell
  } else {
    rowData(d_se)[, "cluster_id"] <- df[, clusterName]
    clusterType <- markersOrCell
  }

  # Calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  rowData(d_counts)[, "cluster_id"] <-
    as.factor(rownames(rowData(d_counts)))

  if (calculateSampleContributionsToClusters) {
    # Transform the cluster cell counts into a plotable format
    counts_df <- assay(d_counts)
    percentage_counts_df <- (counts_df / rowSums(counts_df)) * 100
    t_percentage_counts_df <- t(percentage_counts_df)
    t_percentage_counts_df <- as.data.frame(t_percentage_counts_df)
    colnames(t_percentage_counts_df) <-
      rowData(d_counts)[, "cluster_id"]
    t_percentage_counts_df[, "rownames"] <-
      row.names(t_percentage_counts_df)

    colnames(t_percentage_counts_df) <-
      colnames(t_percentage_counts_df)

    for (columnName in colnames(t_percentage_counts_df)) {
      gc()
      jpeg(file = paste0(
        figureDirectory,
        paste0(
          "sampleContributionTo",
          clusterName,
          columnName,
          "Visits",
          concatinatedVisits,
          clusterType,
          ".jpeg"
        )
      ))
      par(mar = c(1, 1, 1, 1))

      p <-
        ggplot(data = t_percentage_counts_df, aes_string(x = "rownames", y = shQuote(columnName))) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("Patient Sample") + ylab("Percentage") + ggtitle(columnName)
      print(p)
      dev.off()
      gc()

      try(capture.output(
        row.names(t_percentage_counts_df[which(t_percentage_counts_df[, columnName] >=
                                                 samplesContributionToClustersThreshold), ]),
        file = paste0(
          "differentialTestingOutputs/sampleContributionTo",
          clusterName,
          columnName,
          "Visits",
          concatinatedVisits,
          clusterType,
          ".txt"
        )
      ))
    }
  }

  # Calculate cluster medians
  d_medians <- calcMedians(d_se)
  rowData(d_medians)[, "cluster_id"] <-
    as.factor(rownames(rowData(d_counts)))

  # Create experiment design
  design <-
    createDesignMatrix(experimentInfo, cols_design = c("group_id", covariants))

  # Create contrast (the 1 indicates the columns in the design to test)
  contrast <- createContrast(c(0, 1, rep(0, ncol(design) - 2)))

  # Check that design matches control
  nrow(contrast) == ncol(design)
  data.frame(parameters = colnames(design), contrast)

  ### Check this with ethnicity etc
  # Test for differential abundance (DA) of clusters
  res_DA <- testDA_edgeR(d_counts, design, contrast, min_cells = 0)

  # Extract statistics
  res_DA_DT <- as.data.frame(rowData(res_DA))

  # Add - log10(adjusted P-value)
  res_DA_DT[, "minus_log_p_adj"] <- 0 - log10(res_DA_DT[, "p_adj"])
  res_DA_DT[, "diff_expressed"] <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  res_DA_DT$diff_expressed[res_DA_DT$logFC > 0 &
                             res_DA_DT$p_adj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res_DA_DT$diff_expressed[res_DA_DT$logFC < 0 &
                             res_DA_DT$p_adj < 0.05] <- "DOWN"

  differentialAbundanceManhattanPlot(
    res_DA_DT,
    figureDirectory,
    paste0(
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      clusterType,
      "DifferentialAbundanceManhattanPlot.jpeg"
    )
  )

  differentialAbundanceVolcanoPlot(
    res_DA_DT,
    figureDirectory,
    mycolors,
    paste0(
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      clusterType,
      "DifferentialAbundanceVolcanoPlot.jpeg"
    )
  )

  # display table of results for top DA clusters
  topTable(res_DA, format_vals = TRUE)

  # calculate number of significant detected DA clusters at 10% false discovery
  # rate (FDR)
  table(topTable(res_DA, all = TRUE)$p_adj <= differentialAbundanceThreshold)

  # Test for differential states (DS) within clusters
  res_DS <-
    testDS_limma(d_counts,
                 d_medians,
                 design,
                 contrast,
                 plot = FALSE,
                 min_cells = 0)

  # Extract statistics
  res_DS_DT <- as.data.frame(rowData(res_DS))

  # Add - log10(adjusted P-value)
  res_DS_DT[, "minus_log_p_adj"] <- 0 - log10(res_DS_DT[, "p_adj"])
  res_DS_DT[, "id"] <-
    paste0(res_DS_DT[, "cluster_id"], " " , res_DS_DT[, "marker_id"])
  res_DS_DT[, "diff_expressed"] <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  res_DS_DT$diff_expressed[res_DS_DT$logFC > 0 &
                             res_DS_DT$p_adj < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res_DS_DT$diff_expressed[res_DS_DT$logFC < 0 &
                             res_DS_DT$p_adj < 0.05] <- "DOWN"

  differentialStatesVolcanoPlot(
    res_DS_DT,
    figureDirectory,
    mycolors,
    paste0(
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      clusterType,
      "DifferentialStatesVolcanoPlot.jpeg"
    )
  )

  differentialStatesManhattanPlot(
    res_DS_DT,
    figureDirectory,
    paste0(
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      clusterType,
      "DifferentialStatesManhattanPlot.jpeg"
    )
  )

  differentialStatesSaveResults(
    res_DS_DT,
    res_DS,
    res_DA_DT,
    res_DA,
    paste0(
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      clusterType
    )
  )

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })

  gc()
}

performAllDifferentialAbundanceTests <-
  function(directoryName,
           columnNames,
           clusterName,
           markersOrCell) {
    ### Case vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell)

    ### Case vs Controls for first visit for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell)

    ### Fast Progression vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Fast"

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell, progression)

    ### Fast Progression vs Controls for first visit for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell, progression)

    ### Slow Progression vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Slow"

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell, progression)

    ### Slow Progression vs Controls for first visit for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(directoryName, columnNames, clusterName, samplesContributionToClustersThreshold, differentialAbundanceThreshold, calculateSampleContributionsToClusters, group_id, visits, cases, covariants, singleCluster, markersOrCell, progression)

    ### Fast vs Slow for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "fastSlow"
    visits <- c(1)
    cases <- c("Case")
    covariants <-
      c("ageAtVisit",
        "gender",
        "bulbarLimb",
        "ethnicity",
        "timeFromOnsetToVisitInDays")
    singleCluster <- FALSE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Fast vs Slow for first visit for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Bulbar vs Limb for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "bulbarLimb"
    visits <- c(1)
    cases <- c("Case")
    covariants <-
      c("ageAtVisit",
        "gender",
        "fastSlow",
        "ethnicity",
        "timeFromOnsetToVisitInDays")
    singleCluster <- FALSE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Bulbar vs Limb for first visit for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 1 vs Visit 2 for visit 1 & 2 for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "visit"
    visits <- c(1, 2)
    cases <- c("Case")
    covariants <-
      c(
        "ageAtVisit",
        "gender",
        "ethnicity",
        "fastSlow",
        "bulbarLimb",
        "timeFromOnsetToVisitInDays"
      )
    singleCluster <- FALSE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 1 vs Visit 2 for visit 1 & 2 for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 1 vs Visit 3 for visit 1 & 3 for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "visit"
    visits <- c(1, 3)
    cases <- c("Case")
    covariants <-
      c(
        "ageAtVisit",
        "gender",
        "ethnicity",
        "fastSlow",
        "bulbarLimb",
        "timeFromOnsetToVisitInDays"
      )
    singleCluster <- FALSE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 1 vs Visit 3 for visit 1 & 3 for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 2 vs Visit 3 for visit 1 & 3 for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "visit"
    visits <- c(2, 3)
    cases <- c("Case")
    covariants <-
      c(
        "ageAtVisit",
        "gender",
        "ethnicity",
        "fastSlow",
        "bulbarLimb",
        "timeFromOnsetToVisitInDays"
      )
    singleCluster <- FALSE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )


    ### Visit 1 vs Visit 3 for visit 1 & 3 for all cells
    singleCluster <- TRUE

    differentialAbundanceAnalysis(
      directoryName,
      columnNames,
      clusterName,
      samplesContributionToClustersThreshold,
      differentialAbundanceThreshold,
      calculateSampleContributionsToClusters,
      group_id,
      visits,
      cases,
      covariants,
      singleCluster,
      markersOrCell
    )
  }

recalculatePValueAdjustments <-
  function(DA,
           sigCutOff,
           fileNames,
           clusterName,
           markersOrCell,
           flipFoldChange = TRUE) {
    directories <- c("bCells", "monocytes", "tCells", "senescence")

    names(fileNames) <- c("allCells", "cellPopulations")

    for (directory in directories) {
      i <- 1
      for (file in fileNames) {
        names(file) <- names(fileNames)[i]
        i <- i + 1
        filePath <-
          paste0("data/",
                 directory,
                 "/differentialTestingOutputs/",
                 file)
        df <- read.csv(filePath)
        df[, "panel"] <- directory
        if (names(file) == "allCells") {
          if (directory == "bCells") {
            df[, "typeOfCells"] <- "B Cells"
          } else if (directory == "monocytes") {
            df[, "typeOfCells"] <- "Monocytes"
          } else if (directory == "tCells") {
            df[, "typeOfCells"] <- "T Cells"
          } else if (directory == "senescence") {
            df[, "typeOfCells"] <- "Senescent T Cells"
          }
        } else {
          df[, "typeOfCells"] <- df[, "cluster_id"]
        }
        if (exists("combinedDf")) {
          combinedDf <- rbind(combinedDf, df)
        } else {
          combinedDf <- df
        }
      }
    }

    if (flipFoldChange) {
      combinedDf[, "logFC"] <- 0 - combinedDf[, "logFC"]
    }

    # Bonferroni P-Value Adjustment
    combinedDf[, "bonferroni_adjusted_p_val"] <-
      p.adjust(combinedDf[, "p_val"], method = "bonferroni")
    combinedDf[, "minus_log_bonferroni_adjusted_p_val"] <-
      0 - log10(combinedDf[, "bonferroni_adjusted_p_val"])

    # Benjamini and Hochberg P-Value Adjustment
    combinedDf[, "fdr_adjusted_p_val"] <-
      p.adjust(combinedDf[, "p_val"], method = "fdr")
    combinedDf[, "minus_log_fdr_adjusted_p_val"] <-
      0 - log10(combinedDf[, "fdr_adjusted_p_val"])

    # Update differential expression column
    combinedDf$bonferroni_diff_expressed <- "NO"
    combinedDf$bonferroni_diff_expressed[combinedDf[, "logFC"] < 0 &
                                           combinedDf[, "bonferroni_adjusted_p_val"] < sigCutOff] <-
      "DOWN"
    combinedDf$bonferroni_diff_expressed[combinedDf[, "logFC"] > 0 &
                                           combinedDf[, "bonferroni_adjusted_p_val"] < sigCutOff] <-
      "UP"

    combinedDf$fdr_diff_expressed <- "NO"
    combinedDf$fdr_diff_expressed[combinedDf[, "logFC"] < 0 &
                                    combinedDf[, "fdr_adjusted_p_val"] < sigCutOff] <-
      "DOWN"
    combinedDf$fdr_diff_expressed[combinedDf[, "logFC"] > 0 &
                                    combinedDf[, "fdr_adjusted_p_val"] < sigCutOff] <-
      "UP"

    # Update labels
    combinedDf[, "fdr_label"] <- combinedDf[, "typeOfCells"]
    combinedDf$fdr_label[is.na(combinedDf[, "logFC"])] <- NA
    combinedDf$fdr_label[is.na(combinedDf[, "logFC"]) |
                           combinedDf[, "fdr_adjusted_p_val"] > sigCutOff] <-
      NA

    combinedDf[, "bonferroni_label"] <- combinedDf[, "typeOfCells"]
    combinedDf$bonferroni_label[is.na(combinedDf[, "logFC"])] <- NA
    combinedDf$bonferroni_label[is.na(combinedDf[, "logFC"]) |
                                  combinedDf[, "bonferroni_adjusted_p_val"] > sigCutOff] <-
      NA


    # Define Colours
    mycolors <- data.frame(DOWN = "blue",
                           UP = "red",
                           NO = "black")

    setwd("data")

    dir.create("pValueAdjustmentsResults", showWarnings = FALSE)
    dir.create("figures", showWarnings = FALSE)
    figureDirectory <- paste0(getwd(), "/figures/")

    if (DA) {
      jpeg(file = paste0(
        figureDirectory,
        "fdr",
        str_replace_all(str_replace_all(
          str_replace_all(fileNames, " ", ""), ",", ""
        ), "\\.", "")[1],
        markersOrCell,
        ".jpeg"
      ))
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data = combinedDf,
        aes(
          x = logFC,
          y = minus_log_fdr_adjusted_p_val,
          col = fdr_diff_expressed,
          label = fdr_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential Abundance of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()
      gc()

      jpeg(
        file = paste0(
          figureDirectory,
          "bonferroni",
          str_replace_all(str_replace_all(
            str_replace_all(fileNames, " ", ""), ",", ""
          ), "\\.", "")[1],
          markersOrCell,
          ".jpeg"
        )
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data = combinedDf,
        aes(
          x = logFC,
          y = minus_log_bonferroni_adjusted_p_val,
          col = bonferroni_diff_expressed,
          label = bonferroni_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential Abundance of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()
      gc()
    } else {
      # Update marker columsn
      combinedDf[combinedDf[, "marker_id"] == "GPR32...AF488.A", "marker_id"] <-
        "GPR32"
      combinedDf[combinedDf[, "marker_id"] == "GPR32.AF488.A", "marker_id"] <-
        "GPR32"
      combinedDf[combinedDf[, "marker_id"] == "FPRL1...AF647.A", "marker_id"] <-
        "FPRL1"
      combinedDf[combinedDf[, "marker_id"] == "FPRL1.AF647.A", "marker_id"] <-
        "FPRL1"

      # fdr figures
      jpeg(
        file = paste0(
          figureDirectory,
          "gpr32",
          "fdr",
          str_replace_all(str_replace_all(
            str_replace_all(fileNames, " ", ""), ",", ""
          ), "\\.", "")[1],
          markersOrCell,
          ".jpeg"
        )
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data =
          combinedDf[combinedDf[, "marker_id"] == "GPR32", ],
        aes(
          x = logFC,
          y = minus_log_fdr_adjusted_p_val,
          col = fdr_diff_expressed,
          label = fdr_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential States of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()
      gc()

      jpeg(
        file = paste0(
          figureDirectory,
          "fprl1",
          "fdr",
          str_replace_all(str_replace_all(
            str_replace_all(fileNames, " ", ""), ",", ""
          ), "\\.", "")[1],
          markersOrCell,
          ".jpeg"
        )
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data =
          combinedDf[combinedDf[, "marker_id"] == "FPRL1", ],
        aes(
          x = logFC,
          y = minus_log_fdr_adjusted_p_val,
          col = fdr_diff_expressed,
          label = fdr_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential States of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()

      # bonferroni figures
      jpeg(
        file = paste0(
          figureDirectory,
          "gpr32",
          "bonferroni",
          str_replace_all(str_replace_all(
            str_replace_all(fileNames, " ", ""), ",", ""
          ), "\\.", "")[1],
          markersOrCell,
          ".jpeg"
        )
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data =
          combinedDf[combinedDf[, "marker_id"] == "GPR32", ],
        aes(
          x = logFC,
          y = minus_log_bonferroni_adjusted_p_val,
          col = bonferroni_diff_expressed,
          label = bonferroni_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential States of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()
      gc()

      jpeg(
        file = paste0(
          figureDirectory,
          "fprl1",
          "bonferroni",
          str_replace_all(str_replace_all(
            str_replace_all(fileNames, " ", ""), ",", ""
          ), "\\.", "")[1],
          markersOrCell,
          ".jpeg"
        )
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data =
          combinedDf[combinedDf[, "marker_id"] == "FPRL1", ],
        aes(
          x = logFC,
          y = minus_log_bonferroni_adjusted_p_val,
          col = bonferroni_diff_expressed,
          label = bonferroni_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        geom_hline(yintercept = -log10(0.05), col = "red") +
        geom_hline(yintercept = -log10(0.01), col = "red") +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential States of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log Adjusted P-Value")
      print(p)
      dev.off()

      gc()
    }

    write.csv(
      combinedDf,
      paste0(
        "pValueAdjustmentsResults/",
        str_replace_all(str_replace_all(
          str_replace_all(fileNames, " ", ""), ",", ""
        ),
        "\\.",
        "")[1],
        markersOrCell,
        ".csv"
      ),
      row.names = FALSE
    )

    setwd("..")
  }

defineFlowSomCellPopulations <- function(directory, queries) {
  fSOM <-
    readRDS(paste0("data/", directory, "/clusteringOutput/flowSom.rds"))

  cellTypes <- factor(rep("Unlabeled", fSOM$map$nNodes),
                      levels = c("Unlabeled", unique(names(queries))))

  i <- 1

  for (query in queries) {
    query_res <-
      QueryStarPlot(fSOM, query, equalNodeSize = TRUE, plot = FALSE)

    cellTypes[query_res$selected] <- names(queries[i])
    i <- i + 1
  }

  flowsom_cell_populations <- as.factor(fSOM$map$mapping[, 1])
  levels(flowsom_cell_populations) <- cellTypes

  df <-
    read.csv(paste0("data/", directory, '/clusteringOutput/flowSomDf.csv'))
  df[, "flowsom_cell_populations"] <- flowsom_cell_populations

  for (cluster in unique(df[, "clusters_flowsom"])) {
    cell_population <-
      unique(df[df[, "clusters_flowsom"] == cluster, ]$flowsom_cell_populations)
    cell_population <-
      cell_population[cell_population != "Unlabeled"]

    if (length(cell_population) == 1) {
      df[df[, "clusters_flowsom"] == cluster, "flowsom_cell_populations_updated"] <-
        as.character(cell_population)
    } else {
      df[df[, "clusters_flowsom"] == cluster, "flowsom_cell_populations_updated"] <-
        df[df[, "clusters_flowsom"] == cluster, "flowsom_cell_populations"]
    }
  }

  write.csv(df,
            paste0("data/", directory, '/clusteringOutput/flowSomDf.csv'),
            row.names = FALSE)
}

calculateClusterMarkers <-
  function(directoryName,
           clusterName,
           columnNames,
           cutoff,
           markersOrCell = "CellPopulations") {
    if (clusterName == "clusters_flowsom" |
        clusterName == "meta_clusters_flowsom") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/flowSomDf.csv"
        ))
    } else if (clusterName == "clusters_phenograph") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/phenographDf.csv"
        ))
    } else if (clusterName == "clusters_fast_pg") {
      df <-
        read.csv(paste0("data/", directoryName, "/clusteringOutput/fastPGDf.csv"))
    }

    columnNamesMedian <- paste0(columnNames, "_median")
    columnNamesPositive <- paste0(columnNames, "_positive")

    results <-
      data.frame(matrix(
        ncol = 2 * length(columnNamesMedian) + 2,
        nrow = 0
      ))

    for (cluster in unique(df[, clusterName])) {
      new_row <- c(cluster)

      df2 <- df[df[, clusterName] == cluster, ]

      for (column in columnNames) {
        clusterMedian <- median(df2[, column])

        new_row <- append(new_row, clusterMedian)
      }

      new_row <-
        append(new_row, replicate(length(colnames(results)) - length(new_row), NA))

      results <- rbind(new_row, results)
    }

    colnames(results) <-
      c(clusterName,
        columnNamesMedian,
        columnNamesPositive,
        "cell_population")

    i <- 1

    for (column in columnNamesMedian) {
      results[, columnNamesPositive[i]] <-
        results[, columnNamesMedian[i]] > cutoff[i]
      i <- i + 1
    }

    if (markersOrCell == "CellPopulations") {
      cellPopulationMarkers <-
        read.csv(paste0("data/metadata/", directoryName, ".csv"))
    } else if (markersOrCell == "Markers") {
      cellPopulationMarkers <-
        read.csv(paste0("data/metadata/", directoryName, "Markers.csv"))
    }

    cellPopulationMarkers <-
      cellPopulationMarkers[, c("name", columnNamesPositive)]

    if (directoryName == "bCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          IgD...PerCP.Cy5.5.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "IgD...PerCP.Cy5.5.A_positive"] &
            CD24...BV605.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD24...BV605.A_positive"] &
            CD27...BV650.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27...BV650.A_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "monocytes") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          CD11b...17BV421.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b...17BV421.A_positive"] &
            CD11b.activated...PE.Cy7.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b.activated...PE.Cy7.A_positive"] &
            CD14...BV605.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD14...BV605.A_positive"] &
            CD16...PE.CF595.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD16...PE.CF595.A_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "tCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          CD127.BV510.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD127.BV510.A_positive"] &
            CD8.BV650.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8.BV650.A_positive"] &
            CD25.BV786.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD25.BV786.A_positive"] &
            FoxP3.PE.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "FoxP3.PE.A_positive"] &
            CD45RO.PE.CF595.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RO.PE.CF595.A_positive"] &
            CD4.PerCP.Cy5.5.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4.PerCP.Cy5.5.A_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "senescence") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          CD27.BV421.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27.BV421.A_positive"] &
            CD45RA.BV605.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RA.BV605.A_positive"] &
            CD28.BV785.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD28.BV785.A_positive"] &
            KLRG1.PE.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "KLRG1.PE.A_positive"] &
            CD4.PE.CF594.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4.PE.CF594.A_positive"] &
            CCR7.PE.Cy7.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CCR7.PE.Cy7.A_positive"] &
            CD8.PerCP.Cy5.5.A_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8.PerCP.Cy5.5.A_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    }

    write.csv(
      results,
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        markersOrCell,
        ".csv"
      ),
      row.names = FALSE
    )

    write.csv(
      results[is.na(results$cell_population), columnNamesPositive],
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        "Unknown",
        markersOrCell,
        ".csv"
      ),
      row.names = FALSE
    )
  }


calculateVarianceWithinClusters <-
  function(directoryName,
           columnNames,
           clusterName,
           numberOfClusters = NA) {
    if (clusterName == "meta_clusters_flowsom") {
      df <-
        read.csv(
          paste0(
            "data/",
            directoryName,
            "/clusteringOutput/flowSomDf",
            numberOfClusters,
            ".csv"
          )
        )
    } else if (clusterName == "clusters_flowsom") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/flowSomDf.csv"
        ))
    } else if (clusterName == "clusters_phenograph") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/phenographDf.csv"
        ))
    } else if (clusterName == "clusters_fast_pg") {
      df <-
        read.csv(paste0("data/", directoryName, "/clusteringOutput/fastPGDf.csv"))
    }

    df <- df[, c(columnNames, clusterName)]

    #Function to calculate the total within-cluster variance for a model of k-classes
    #df is a dataframe expected to contain only the columns used for clustering and a column indicating the allocated cluster (titled "clusters_flowsom")

    calcVar <- function(df2) {
      #Compute total variance across all clusters and predictor variables
      totvar <- df2 %>%
        group_by_at(clusterName) %>%                                 #Split clusters
        dplyr::summarise(across(everything(), ~ var(.))) %>%      #Identify variance for each variable by cluster
        ungroup() %>%                                             #Remove grouping by cluster
        select(-paste(clusterName)) %>%                             #drop cluster identifying column
        sum()                                                     #sum all variances

      return(totvar)

    }

    frame <- list()
    for (cluster in unique(df[, clusterName])) {
      clusterDf <- df[df[, clusterName] == cluster, ]
      frame <- append(frame, list(clusterDf))
    }

    #Apply calcVar across a list of data.frames which individually detail each cluster allocation
    clustVar <- sapply(frame, calcVar)

    saveRDS(
      clustVar,
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        numberOfClusters,
        "flowSomDf.rds"
      )
    )

    figureDirectory <- paste0("data/", directoryName, "/figures/")

    jpeg(file = paste0(
      figureDirectory,
      "elbowPlot",
      clusterName,
      numberOfClusters,
      ".jpeg"
    ))
    #Plot the elbow plot across variances from clustVar
    elbow_plot <-
      ggplot(data.frame(clustVar),
             aes(
               x = as.factor(1:length(clustVar)),
               y = clustVar,
               group = 1
             )) +
      geom_point() +
      geom_line() +
      xlab("Number of clusters") +
      ylab("Sum of within-cluster variance") +
      theme_bw()
    try(print(elbow_plot))
    dev.off()
  }

consolidateFlowSomClusters <-
  function(directoryName,
           columnNames,
           clusterName,
           numberOfClusters) {
    df <-
      read.csv(
        paste0(
          "data/",
          directoryName,
          "/clusteringOutput/",
          "flowSomDf",
          numberOfClusters[1],
          ".csv"
        )
      )

    df <- df[, c(columnNames)]

    for (number in numberOfClusters) {
      x <-
        read.csv(
          paste0(
            "data/",
            directoryName,
            "/clusteringOutput/",
            "flowSomDf",
            number,
            ".csv"
          )
        )

      x <- list(x[, clusterName])
      names(x) <- paste0(clusterName, number)

      df <- cbind(df, x)
    }

    write.csv(
      df,
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        "flowSomDfAllMetaClusters.csv"
      ),
      row.names = FALSE
    )
  }

elbowPlot <-
  function(directoryName,
           columnNames,
           numberOfClusters) {
    df <-
      read.csv(
        paste0(
          "data/",
          directoryName,
          "/clusteringOutput/",
          "flowSomDfAllMetaClusters.csv"
        )
        #, nrows = 100000
      )

    head(df)

    data <- df[, columnNames]
    data[, "ID"] <- row.names(data)
    data <- data[, c("ID", columnNames)]
    head(data)

    clust <- df[,!colnames(df) %in% columnNames]
    colnames(clust) <- numberOfClusters
    clustColnames <- colnames(clust)
    clust[, "ID"] <- row.names(clust)
    clust[, "1"] <- 1
    clust <- clust[, c("ID", "1", clustColnames)]
    head(clust)

    calcVar <- function(clust, data, clust_ID) {
      clust_wIDs <-
        data.frame(clust_ID, clust) #Recombine the clust_ID and clust vectors into a dataframe
      colnames(clust_wIDs) <-
        c("ID", "clust") #Rename to appropriate names

      colnames(data)[1] <-
        "ID"               #Rename first column of main dataset to be "ID", since this will be matched within left-join


      #Join the cluster assignments to the main dataset by the ID index
      df <- data %>%
        left_join(clust_wIDs, by = "ID")

      #Compute total variance across all clusters and predictor variables
      totvar <- df %>%
        select(-ID) %>%                                           #Drop the ID column
        group_by(clust) %>%                                       #group by clusters
        dplyr::summarise(across(everything(), ~ var(.))) %>%      #Identify variance for each variable by cluster
        ungroup() %>%                                             #Remove grouping by cluster
        select(-clust)                                       #drop cluster identifying column

      totvar <- sum(totvar ^ 2) #square and sum all variance

      return(totvar)                                              #Return total within-cluster variances for this clustering solution

    }

    #clust_ID <- clust[1]
    #clust <- clust[2]

    clustVar <-
      sapply(clust[-1], calcVar, data = data, clust_ID = clust[1])

    #Generate factor variable indicating the x-axis labels from column names from clustering data
    #Factor ensures correct ordering
    xlevs <-
      factor(colnames(clust[-1]), levels = colnames(clust[-1]))

    #Plot the elbow plot across variances from clustVar
    elbow_plot <-
      ggplot(data.frame(clustVar), aes(x = xlevs, y = clustVar, group = 1)) +
      geom_point() +
      geom_line() +
      xlab("Number of clusters") +
      ylab("Sum of squared within-cluster variance") +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5,
        hjust = 1
      )) #Orient x axis label to allow for longer names, hjust=1 is right aligned, vjust is center aligned

    ggsave(
      plot = elbow_plot,
      filename = paste0("data/", directoryName, "/figures/Elbow_plot.pdf"),
      units = "mm",
      width = 200,
      height = 150,
      device = cairo_pdf
    )
  }


generateHeatmap <-
  function(clusterName, directoryName, columnNames, markersOrCell, markerType, prettyColumnNames) {
    # clusterName <- "meta_clusters_flowsom"
    # directoryName <- "bCells"
    # columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
    # markersOrCell <- "Markers"
    if (clusterName == "meta_clusters_flowsom") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/flowSomDf.csv"
        ))
    } else if (clusterName == "clusters_flowsom") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/flowSomDf.csv"
        ))
    } else if (clusterName == "clusters_phenograph") {
      df <-
        read.csv(paste0(
          "data/",
          directoryName,
          "/clusteringOutput/phenographDf.csv"
        ))
    } else if (clusterName == "clusters_fast_pg") {
      df <-
        read.csv(paste0("data/", directoryName, "/clusteringOutput/fastPGDf.csv"))
    }

    cellPopulationMarkers <-
      read.csv(
        paste0(
          "data/",
          directoryName,
          "/clusteringOutput/",
          clusterName,
          "CellPopulations.csv"
        )
      )

    markerPopulationMarkers <-
      read.csv(paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        "Markers.csv"
      ))

    colnames(markerPopulationMarkers)[ncol(markerPopulationMarkers)] <-
      "marker_population"

    populations <-
      merge(cellPopulationMarkers, markerPopulationMarkers, by = clusterName)

    populations <-
      populations[, c(clusterName, "cell_population", "marker_population")]

    df <- df[, c(columnNames, clusterName)]

    colnames(df)[1:length(columnNames)] <- prettyColumnNames

    colnames(df)[ncol(df)] <- "clusters"

    figureDirectory <- paste0("data/", directoryName, "/figures/")

    df <- df[df[, "clusters"] %in% populations[, clusterName],]
    populations <- populations[populations[, clusterName] %in% df[, "clusters"],]


    jpeg(file = paste0(figureDirectory,
                       "na.jpeg"))
    metacluster_colours <-
      as.vector(colortools::wheel('#ce302e', num = nrow(populations)))
    dev.off()

    names(metacluster_colours) <- populations[, if (markersOrCell == "CellPopulations") {
      "cell_population"
    } else {
      "marker_population"
    }]

    summarisedDf <- df %>%
      group_by(clusters) %>%
      summarise_all(median) %>%
      remove_rownames() %>%
      column_to_rownames('clusters') %>%
      as.matrix() %>%
      pheatmap:::scale_rows()

    row_order <-
      populations[match(rownames(summarisedDf), populations[, clusterName]),
                  if (markersOrCell == "CellPopulations") {
                    "cell_population"
                  } else {
                    "marker_population"
                  }]

    left_anno <-
      HeatmapAnnotation(
        Celltype = row_order,
        which = 'row',
        col = list('Celltype' = metacluster_colours),
        show_legend = F,
        show_annotation_name = F
      )

    fig <-
      Heatmap(
        summarisedDf,
        rect_gp = gpar(col = "white", lwd = 2),
        column_names_rot = 45,
        show_heatmap_legend = T,
        row_names_side = 'left',
        heatmap_legend_param = list(title = 'Z-score'),
        left_annotation = left_anno,
        row_split = row_order,
        row_gap = unit(2, "mm"),
        border = T,
        row_title_rot = 0,
        col = circlize::colorRamp2(seq(-4, 4, length = 3), c('#4575B4', 'white', '#D73027'))
      )

    png(
      paste0(
        figureDirectory,
        clusterName,
        markersOrCell,
        markerType,
        "heatmap.jpeg"
      ),
      width = 4000,
      height = 2000,
      res = 320
    )

    print(fig)
    dev.off()
  }
