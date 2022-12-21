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
  install.packages(c("doSNOW", "doParallel", "doMPI", "pheatmap"))

  gc()
}

loadlibraries <- function() {
  try(library(pkgconfig))
  try(library(data.table))
  try(library(statmod))
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
  try(library(cytofkit2))
  try(library(ggrepel))
  try(library(tidyverse))
  try(library(diffcyt))
  try(library(corrplot))
  try(library(cluster))
  try(library(factoextra))
  try(library(readxl))
  try(library(fda))
}

ungzipFiles <- function() {
  workingDirectory <- getwd()

  dataDirectorys <- c("/data/gpr32BCells",
                      "/data/gpr32BMonocytes",
                      "/data/gpr32BSenescence",
                      "/data/gpr32TCells")

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

  dataDirectorys <- c("/data/gpr32BCells",
                      "/data/gpr32BMonocytes",
                      "/data/gpr32BSenescence",
                      "/data/gpr32TCells")

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

rawDensityPlots <- function(directoryName,
                          columnNames,
                          prettyColumnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

    # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  dfs <- lapply(filenames, read.csv,
                   sep = ",",
                   header = TRUE,
                   stringsAsFactors = FALSE)

  df <- bind_rows(dfs)

  df <- df[, columnNames]

  colnames(df) <- prettyColumnNames

  dir.create("figures", showWarnings = FALSE)

  figureDirectory <- paste0(getwd(), "/figures/")

  for (columnName in prettyColumnNames) {
    jpeg(file = paste0(figureDirectory,
                       "rawCombinedDensityPlot",
                       columnName,
                       ".jpeg"), quality = 100)
    d <- density(df[, columnName])
    print(plot(d))
    dev.off()
  }

  i <- 1
  for (columnName in prettyColumnNames) {
    jpeg(file = paste0(figureDirectory,
                       "linearPlusLogCombinedDensityPlot",
                       columnName,
                       ".jpeg"), quality = 100)
    x <- df[, columnName]
    x[x > automatedcofactors[i]] <- asinh(x[x > automatedcofactors[i]]) + automatedcofactors[i]

    d <- density(x)
    print(plot(d))
    dev.off()
  }

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

preprocessing <- function(directoryName,
                          columnNames,
                          prettyColumnNames,
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

  filenames <- filenames[!filenames %in% c("BLT00254_T2.csv",
                                 "BLT00271_T2.csv",
                                 "BLT00275_T2.csv",
                                 "BLT00286_T2.csv",
                                 "HC_BLT00282_T2.csv",
                                 "HC_BLT00287_T2.csv")]

  if (test) {
    filenames <- filenames[1:4]
  }

  read.flow_csv <- function(pathIN) {
    message(pathIN)
    raw <-
      fread(
        pathIN#,
        # ,
        # sep = ",",
        # header = TRUE,
        # stringsAsFactors = FALSE,
         #nrows = 100
      )
    raw <- as.data.frame(raw)
    # message(colnames(raw))
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

    #if ("GPR18...AF488.A" %in% colnames(raw)){
    #  print("")
    #  print(pathIN)
    #  #print(colnames(raw))
    #}
    return(unique(raw))
  }

  # Read all:
  dfs <- sapply(filenames, read.flow_csv, simplify = FALSE)

  xlist <- c()

  for (x in seq(length(dfs))) {
    df <- dfs[[x]]
    # print("")
    # print(colnames(df))
    # print("")

    xlist <- append(xlist, colnames(df))
    xcolumnNames <- columnNames[columnNames %in% colnames(df)]
    xprettyColumnNames <- prettyColumnNames[columnNames %in% colnames(df)]

    # print(columnNames[!columnNames %in% colnames(df)])

    df <- df[, xcolumnNames]
    colnames(df) <- xprettyColumnNames
    dfs[[x]]  <- df

    if(length(xcolumnNames) < 8){
      print(names(dfs[x]))
      print(xcolumnNames)
    }
  }

  unique(xlist)

  ##############################
  #REWRITE TO FLOWFRAME/FLOWSET#
  ##############################

  ## Defining a function to rewrite a csv into a flowframe:
  csv_2_ff <- function(dat) {
    message(dimnames(dat)[[2]])
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
  dfs_ff <- sapply(dfs, function(x)
    csv_2_ff(x), simplify = FALSE)
  gc()
  #rm(dfs)

  # rewrite to flowset
  dfs_fs <- as(dfs_ff, "flowSet")
  gc()
  # rm(dfs_ff)

  dir.create("figures", showWarnings = FALSE)

  prettyColumnNames <- unique(prettyColumnNames)

  if (is.null(automatedcofactors)) {
    gc()
    jpeg(file = "figures/automatedcofactors.jpeg", quality = 100)
    automatedcofactors <- estParamFlowVS(dfs_fs, prettyColumnNames)
    dev.off()
    try(capture.output(automatedcofactors,
                       file = "dataPPOutput/automatedcofactors.txt"))

    gc()
  }

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels = prettyColumnNames,
                               cofactor = automatedcofactors)
  gc()
  # rm(automatedcofactors)


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
    warpSet(dfs_fs_t_auto, stains = prettyColumnNames)
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

  # columnName <- prettyColumnNames[1]

  for (columnName in prettyColumnNames) {
    gc()
    columnNameFormula <- as.formula(paste(" ~ ", columnName))
    gc()

    gc()
    jpeg(file = paste0(
      figureDirectory,
      "rawDensityPlot",
      str_replace_all(columnName, "\\.", ""),
      ".jpeg"
    ), quality = 100)
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
    ), quality = 100)
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
      ), quality = 100
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
      ), quality = 100
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

convertToDataFrame <- function(directoryName, columnNames) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  figureDirectory <- paste0(getwd(), "/figures/")

  dirFCS <- paste0(getwd(), "/dataPPOutput")

  pathST <- "X:/Users/guypw/OneDrive/Documents.txt"

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

  fwrite(df, 'dataPPOutput/rawDf.csv', row.names = FALSE)
  gc()

  updatedColumnNames <- append(columnNames, "fileName")

  df <- df[, updatedColumnNames]
  gc()
  fwrite(df, 'dataPPOutput/columnsOfInterestDf.csv', row.names = FALSE)
  gc()

  for (col in columnNames) {
    jpeg(file = paste0(figureDirectory,
                       "rawDensityPlotCol",
                       col,
                       ".jpeg"), quality = 100)
    d <- density(df[, col])
    try(print(plot(d, main = col)))
    dev.off()
  }

  # Scale data
  scaledDf <- df

  for (col in columnNames) {
    ninetyNinthQuantile <-
      quantile(scaledDf[, col], probs = c(0.001, 0.999))
    scaledDf <-
      scaledDf[scaledDf[, col] >= min(ninetyNinthQuantile) &
                 scaledDf[, col] <= max(ninetyNinthQuantile),]
  }

  scaledDf[, columnNames[!columnNames %in% c("GPR32", "GPR18")]] <- scale(scaledDf[, columnNames[!columnNames %in% c("GPR32", "GPR18")]])
  fwrite(scaledDf, 'dataPPOutput/scaledDf.csv', row.names = FALSE)

  for (col in columnNames) {
    jpeg(file = paste0(figureDirectory,
                       "scaledDensityPlotCol",
                       col,
                       ".jpeg"), quality = 100)
    d <- density(scaledDf[, col])
    try(print(plot(d, main = col)))
    dev.off()
  }

  minMaxDf <- df

  for (col in columnNames) {
    dfMax <- max(minMaxDf[, col])
    dfMin <- min(minMaxDf[, col])

    scaledColumn <- (minMaxDf[, col] - dfMin) / (dfMax - dfMin)

    minMaxDf[, col] <- scaledColumn

    jpeg(file = paste0(figureDirectory,
                       "minMaxDensityPlotCol",
                       col,
                       ".jpeg"), quality = 100)
    d <- density(scaledColumn)
    print(plot(d))
    dev.off()
  }

  fwrite(minMaxDf, 'dataPPOutput/minMaxScaledDf.csv', row.names = FALSE)

  outliersMinMaxDf <- df

  # remove outliers
  for (col in columnNames) {
    ninetyNinthQuantile <-
      quantile(outliersMinMaxDf[, col], probs = c(0.001, 0.999))
    outliersMinMaxDf <-
      outliersMinMaxDf[outliersMinMaxDf[, col] >= min(ninetyNinthQuantile) &
                         outliersMinMaxDf[, col] <= max(ninetyNinthQuantile),]
  }

  columnNames <- columnNames[!columnNames %in% c("GPR32", "GPR18")]

  # plot out puts
  for (col in columnNames) {
    dfMax <- max(outliersMinMaxDf[, col])
    dfMin <- min(outliersMinMaxDf[, col])

    scaledColumn <-
      (outliersMinMaxDf[, col] - dfMin) / (dfMax - dfMin)

    outliersMinMaxDf[, col] <- scaledColumn

    jpeg(file = paste0(
      figureDirectory,
      "outliersRemovedMinMaxDensityPlotCol",
      col,
      ".jpeg"
    ), quality = 100)
    d <- density(scaledColumn)
    print(plot(d))
    dev.off()
  }

  fwrite(outliersMinMaxDf,
            'dataPPOutput/outliersRemoveMinMaxScaledDf.csv',
            row.names = FALSE)


  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

convertToFCS <- function(directoryName) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName, "/dataPPOutput"))

  df <- fread(file='outliersRemoveMinMaxScaledDf.csv')
  df <- as.data.frame(df)

  dfs <- list()

  fileNames <- unique(df$fileName)

  for (fileName in fileNames) {
    x <- df[df$fileName == fileName, ]
    x <- x[,!(colnames(x) == "fileName")]
    x <- list(x)
    names(x) <- fileName

    dfs <- append(dfs, x)
  }

  csv_2_ff <- function(dat) {
    # Compute required metadata - column names with description -
    # ranges, min, and max settings
    meta <- data.frame(
      name = names(dat),
      desc = names(dat),
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

  dfs_fs <- as(dfs_ff, "flowSet")
  gc()

  dir.create("scaledFcs", showWarnings = FALSE)

  write.flowSet(dfs_fs,
                outdir = 'scaledFcs',
                filename = paste0(sampleNames(dfs_fs), ".fcs"))

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
    setwd("..")
  })
}

copyToClusteringOutput <- function(directoryName) {
  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  df <- fread(file="./dataPPOutput/scaledDf.csv")
  df <- as.data.frame(df)

  dir.create("clusteringOutput", showWarnings = FALSE)

  fwrite(df, 'clusteringOutput/clusteringOutputs.csv', row.names = FALSE)

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}

read.twoflowdat <- function(dir, path_CSPLR_ST = "") {
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


flowsomClustering <-
  function(directoryName,
           columnNames,
           numberOfClusters,
           df,
           fcsDf) {
    workingDirectory <- getwd()

    setwd(paste0("./data/", directoryName))

    dirFCS <- paste0(getwd(), "/dataPPOutput/scaledFcs")

    columnIndexes <- seq(length(columnNames))

    seed <- 100

    #run flowsom
    flowsom <- FlowSOM(
      input = dirFCS,
      transform = FALSE,
      scale = FALSE,
      colsToUse = columnIndexes,
      nClus = numberOfClusters,
      seed = seed
    )

    # Get metaclustering per cell
    clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
    meta_clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
    levels(meta_clusters_flowsom) <- flowsom$metaclustering

    #add flowsom clusters to dataframe
    df[, "clusters_flowsom"] <- clusters_flowsom
    df[, paste0("meta_clusters_flowsom", as.character(numberOfClusters))] <-
      meta_clusters_flowsom

    fwrite(df,
              'clusteringOutput/clusteringOutputs.csv',
              row.names = FALSE)
    try(FlowSOMmary(flowsom,
                    plotFile = paste0("clusteringOutput/FlowSOMmary", numberOfClusters, ".pdf")))

    tryCatch({
      setwd(workingDirectory)
    },
    error = function(cond) {
      setwd("..")
      setwd("..")
    })

    return(df)
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

  df <- fread(file='clusteringOutput/clusteringOutputs.csv')
  df <- as.data.frame(df)

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
  df[, "clusters_phenograph"] <- clusters_phenograph

  fwrite(df, 'clusteringOutput/clusteringOutputs.csv', row.names = FALSE)

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

  df <- fread(file='clusteringOutput/clusteringOutputs.csv')
  df <- as.data.frame(df)
  #df <- fread(file='clusteringOutput/flowsomClusterStability.csv')

  gc()

  fastPGResults <-
    FastPG::fastCluster(as.matrix(df[, columnNames]), knn, 5)

  clusters_fast_pg <- fastPGResults$communities

  #add phenograph clusters to expression data frame
  df[, "clusters_fast_pg"] <- clusters_fast_pg

  fwrite(df, 'clusteringOutput/clusteringOutputs.csv', row.names = FALSE)

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

  df <-
    fread(file='clusteringOutput/clusteringOutputs.csv')
  df <- as.data.frame(df)

  try({randomNumbers <-
    sample(seq(nrow(df)), 50000, replace = FALSE)

  df <- df[randomNumbers,]})

  umap <- umap(df[, columnNames],
               n_neighbors = knn,
               min_dist = 0.001,
               verbose = TRUE)
  umap <- as.data.frame(umap)
  colnames(umap) <-
    c('umap_1', 'umap_2')
  df <- cbind(df, umap)

  fwrite(df, 'clusteringOutput/umapDf.csv', row.names = FALSE)
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

  df <- fread(file='clusteringOutput/umapDf.csv')
  df <- as.data.frame(df)

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
    ), quality = 100)
    plot <- viz.umap(dat = df, param.name = columnName)
    try(print(plot))
    dev.off()
    gc()
  }

  metaFlowSomMarkers <-
    fread(file='clusteringOutput/meta_clusters_flowsomMarkers.csv')
  metaFlowSomMarkers <- as.data.frame(metaFlowSomMarkers)

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
    fread(file='clusteringOutput/clusters_flowsomMarkers.csv')
  flowSomMarkers <- as.data.frame(flowSomMarkers)
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

  try({
    fastPgMarkers <-
      fread(file='clusteringOutput/clusters_fast_pgMarkers.csv')
    fastPgMarkers <- as.data.frame(fastPgMarkers)
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
  })

  try({
    phenographMarkers <-
      fread(file='clusteringOutput/clusters_phenographMarkers.csv')
    phenographMarkers <- as.data.frame(phenographMarkers)
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
  })



  metaFlowSomcellPopulations <-
    fread(file='clusteringOutput/meta_clusters_flowsomCellPopulations.csv')
  metaFlowSomcellPopulations <- as.data.frame(metaFlowSomcellPopulations)
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

  try({
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
  })

  try({
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
  })

  #visualize and label clusters on umap
  gc()
  label_flowsom_umap <- df %>% group_by(flowsom_cell_population) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_flowsom_cell_population) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)

  try({
    label_pheno_umap <-
      df %>% group_by(phenograph_cell_population) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })

  try({
    label_fastpg_umap <- df %>% group_by(fastpg_cell_population) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomCellPopulations.jpeg"), quality = 100)
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
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomCellPopulations.jpeg"), quality = 100)
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

  try({
    jpeg(file = paste0(figureDirectory, "umapPhenographCellPopulations.jpeg"), quality = 100)
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
  })
  dev.off()

  gc()

  try({
    jpeg(file = paste0(figureDirectory, "umapFastPGCellPopulations.jpeg"), quality = 100)
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
  })
  dev.off()

  label_flowsom_umap <- df %>% group_by(clusters_flowsom) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_clusters_flowsom) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  try({
    label_pheno_umap <- df %>% group_by(clusters_phenograph) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })
  try({
    label_fastpg_umap <- df %>% group_by(clusters_fast_pg) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomClusters.jpeg"), quality = 100)
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
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomClusters.jpeg"), quality = 100)
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
  try({
    jpeg(file = paste0(figureDirectory, "umapPhenographClusters.jpeg"), quality = 100)
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
  })
  dev.off()
  gc()
  try({
    jpeg(file = paste0(figureDirectory, "umapFastPGClusters.jpeg"), quality = 100)
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
  })
  dev.off()

  label_flowsom_umap <- df %>% group_by(flowsom_markers) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  label_meta_flowsom_umap <-
    df %>% group_by(meta_flowsom_markers) %>%
    dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  try({
    label_pheno_umap <- df %>% group_by(phenograph_markers) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })
  try({
    label_fastpg_umap <- df %>% group_by(fastpg_markers) %>%
      dplyr::select(umap_1, umap_2) %>% summarize_all(mean)
  })

  gc()
  jpeg(file = paste0(figureDirectory, "umapFlowsomMarkers.jpeg"), quality = 100)
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
  jpeg(file = paste0(figureDirectory, "umapMetaFlowsomMarkers.jpeg"), quality = 100)
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
  try({
    jpeg(file = paste0(figureDirectory, "umapPhenographMarkers.jpeg"), quality = 100)
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
  })
  dev.off()
  gc()
  try({
    jpeg(file = paste0(figureDirectory, "umapFastPGMarkers.jpeg"), quality = 100)
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
  })
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

    df <-
      fread(file='clusteringOutput/clusteringOutputs.csv')
    df <- as.data.frame(df)

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

    fwrite(df, 'clusteringOutput/diffusionMapDf.csv', row.names = FALSE)
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
    ), quality = 100)
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
  try({
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
  })

  try({
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
  })


  #visualize and label clusters on umap
  gc()
  label_flowsom_dm <- df %>% group_by(flowsom_cell_population) %>%
    dplyr::select(DC1, DC2) %>% summarize_all(mean)
  label_meta_flowsom_dm <-
    df %>% group_by(meta_flowsom_cell_population) %>%
    dplyr::select(DC1, DC2) %>% summarize_all(mean)
  try({
    label_pheno_dm <- df %>% group_by(phenograph_cell_population) %>%
      dplyr::select(DC1, DC2) %>% summarize_all(mean)
  })
  try({
    label_fastpg_dm <- df %>% group_by(fastpg_cell_population) %>%
      dplyr::select(DC1, DC2) %>% summarize_all(mean)
  })

  gc()
  jpeg(file = paste0(
    figureDirectory,
    "diffusionMapFlowsomCellPopulations.jpeg"
  ), quality = 100)
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
  ), quality = 100)
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
  try({
    jpeg(file = paste0(
      figureDirectory,
      "diffusionMapPhenographCellPopulations.jpeg"
    ), quality = 100)
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
  })
  dev.off()
  gc()
  try({
    jpeg(file = paste0(
      figureDirectory,
      "diffusionMapFastPGCellPopulations.jpeg"
    ), quality = 100)
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
  })
  dev.off()

  label_flowsom_dm <- df %>% group_by(clusters_flowsom) %>%
    dplyr::select(DC1, DC2) %>% summarize_all(mean)
  label_meta_flowsom_dm <-
    df %>% group_by(meta_clusters_flowsom) %>%
    dplyr::select(DC1, DC2) %>% summarize_all(mean)
  try({
    label_pheno_dm <- df %>% group_by(clusters_phenograph) %>%
      dplyr::select(DC1, DC2) %>% summarize_all(mean)
  })
  try({
    label_fastpg_dm <- df %>% group_by(clusters_fast_pg) %>%
      dplyr::select(DC1, DC2) %>% summarize_all(mean)
  })

  gc()
  jpeg(file = paste0(figureDirectory, "diffusionMapFlowsomClusters.jpeg"), quality = 100)
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
  jpeg(file = paste0(figureDirectory, "diffusionMapMetaFlowsomClusters.jpeg"), quality = 100)
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
  try({
    jpeg(file = paste0(figureDirectory, "diffusionMapPhenographClusters.jpeg"), quality = 100)
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
  })
  dev.off()
  gc()
  try({
    jpeg(file = paste0(figureDirectory, "diffusionMapFastPGClusters.jpeg"), quality = 100)
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
  })
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
    ), quality = 100)
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

  jpeg(file = paste0(figureDirectory, "mclustBic.jpeg"), quality = 100)
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
  jpeg(file = paste0(figureDirectory, "apcluster1.jpeg"), quality = 100)
  par(mar = c(1, 1, 1, 1))
  p <- heatmap(d.apclus)
  print(p)
  dev.off()
  gc()

  jpeg(file = paste0(figureDirectory, "apcluster1.jpeg"), quality = 100)
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
                     figureTitle), quality = 100)
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
                     figureTitle), quality = 100)
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
                     figureTitle), quality = 100)
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
                     figureTitle), quality = 100)
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
  fwrite(
    res_DS_DT,
    paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialStatesStatistics.csv'
    ),
    row.names = FALSE
  )
  fwrite(
    res_DA_DT,
    paste0(
      'differentialTestingOutputs/',
      datasetTitle,
      'DifferentialAbundanceStatistics.csv'
    ),
    row.names = FALSE
  )
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

  fwrite(df, 'dataPPOutput/gatedDf.csv', row.names = FALSE)

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

    fwrite(df, 'dataPPOutput/gatedDf.csv', row.names = FALSE)

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

differentialAbundanceAnalysis <- function(df,
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
                                          markersOrCell,
                                          progression = "",
                                          siteOfOnset = "",
                                          blocking = NULL) {
  concatinatedVisits <- toString(visits)

  # Read experiment data
  experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

  experimentInfo <- as.data.frame(experimentInfo)

  experimentInfo <- experimentInfo[experimentInfo$experiment == "flowCytometry", ]

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

  experimentInfo <- updateClinicalData(experimentInfo, directoryName)

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

  if (siteOfOnset != "") {
    experimentInfo <-
      experimentInfo[experimentInfo[, "BulbarLimb"] == siteOfOnset |
                       experimentInfo[, "caseControl"] == "Control", ]
  }
  experimentInfo[, "group_id"] <- NA
  experimentInfo <-
    experimentInfo[!is.na(experimentInfo[, "ethnicity"]),]
  experimentInfo[, "ageAtVisitDouble"] <-
    as.double(experimentInfo[, "ageAtVisit"])
  experimentInfo[, "timeFromVisit1InYearsDouble"] <-
    as.double(experimentInfo[, "timeFromVisit1InYears"])
  experimentInfo[, "visit"] <-
    as.double(experimentInfo[, "visit"])

  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")

  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dir.create("differentialTestingOutputs", showWarnings = FALSE)
  figureDirectory <- paste0(getwd(), "/figures/")

  df <- df[order(df[, "fileName"]), ]

  # Filter df for samples
  df <- df[df[, "fileName"] %in% experimentInfo[, "sample_id"], ]

  colnames(experimentInfo)[1] <- "fileID"

  if (markersOrCell != "Clusters") {
    df <-
      merge(df,
            cellPopulationMarkers[, c(clusterName, "cell_population")],
            by = clusterName,
            all.x = TRUE)

  }

  df <-
    merge(df,
          experimentInfo[, c("sample_id", "fileID")],
          by.x = "fileName",
          by.y = "sample_id",
          all.x = TRUE)

  # Filter df for samples
  experimentInfo <-
    experimentInfo[experimentInfo[, "sample_id"] %in% df[, "fileName"],]

  # Reorder df and experiemnt info
  df <- df[order(df[, "fileName"]), ]
  experimentInfo <-
    experimentInfo[order(experimentInfo[, "sample_id"]), ]

  experimentInfo[, "group_id"] <-
    factor(experimentInfo[, group_id])
  experimentInfo[, "patient_id"] <-
    factor(experimentInfo[, "patient_id"])
  experimentInfo[, "fastSlow"] <-
    factor(experimentInfo[, "fastSlow"])
  experimentInfo[, "BulbarLimb"] <-
    factor(experimentInfo[, "BulbarLimb"])
  experimentInfo[, "caseControl"] <-
    factor(experimentInfo[, "caseControl"])
  experimentInfo[, "ethnicity"] <-
    factor(experimentInfo[, "ethnicity"])
  experimentInfo[, "patient_id"] <-
    factor(experimentInfo[, "patient_id"])
  experimentInfo[, "sample_id"] <-
    factor(experimentInfo[, "sample_id"])
  experimentInfo[, "gender"] <-
    factor(experimentInfo[, "gender"])

  # Check df and experiment data are in the same order
  message(all(unique(df[, "fileName"]) == unique(experimentInfo[, "sample_id"])))

  # Extract only the relevant columns
  minimalDf <- df[, c(clusterName, "fileID", columnNames)]

  # split the dataframe into a dataframe for each file
  listOfDfs <- list()
  for (file in unique(minimalDf[, "fileName"])) {
    minimalDfExtract <- minimalDf[minimalDf[, "fileName"] == file, ]
    minimalDfExtract <-
      minimalDfExtract[,!(names(minimalDfExtract) %in% c("fileName"))]
    listOfDfs <- append(listOfDfs, list(minimalDfExtract))
  }

  # Create marker information
  markerColumnNames <- c(clusterName, "fileID", columnNames[columnNames != "fileName"])
  markerInformation <- data.frame(markerColumnNames)
  markerInformation[, "channel_name"] <- markerColumnNames
  markerInformation[, "marker_name"] <- markerColumnNames
  markerInformation[, "marker_class"] <-
    c(rep("none", 2),
      rep("type", length(columnNames) - 2),
      rep("state", 1))
  markerInformation <-
    markerInformation[, 2:ncol(markerInformation)]
  markerInformation[, "marker_class"] <- as.factor(
    markerInformation[, "marker_class"])

  # Transform the input into the correct format
  d_se <- prepareData(listOfDfs, experimentInfo, markerInformation)
  rowData(d_se)[, "file_ID"] <- assay(d_se)[, "fileID"]
  rowData(d_se)[, "row_ID"] <- rownames(rowData(d_se))

  head(rowData(d_se))

  # Test Experiment data was created successfully
  message(all(rowData(d_se)[, "file_ID"] == assay(d_se)[, "fileID"]))

  if (singleCluster) {
    rowData(d_se)[, "cluster_id"] <- 1
    clusterType <- "AllCells"
  } else if (markersOrCell != "Clusters") {
    rowData(d_se)[, "clusterID"] <- assay(d_se)[, clusterName]
    rowData(d_se) <- merge(rowData(d_se),
                           cellPopulationMarkers[, c(clusterName, "cell_population")],
                           by.x = "clusterID", by.y = clusterName,
                           all.x = TRUE)
    rowData(d_se)[, "cluster_id"] <- rowData(d_se)[, "cell_population"]

    #all(!rowData(d_se)[, "cell_population"] %in% cutoffDf[, clusterName])

    rowData(d_se) <- rowData(d_se)[order(as.integer(
                                         rowData(d_se)[, "row_ID"])), ]

    clusterType <- markersOrCell
  } else {
    rowData(d_se)[, "cluster_id"] <- assay(d_se)[, clusterName]
    clusterType <- markersOrCell
  }

  message(all(rowData(d_se)[, "file_ID"] == assay(d_se)[, "fileID"]))

  rowData(d_se)[, "cluster_id"] <- as.factor(rowData(d_se)[, "cluster_id"])

  # Calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  #rowData(d_counts)[, "cluster_id"] <-
  #  as.factor(rownames(rowData(d_counts)))

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
      ), quality = 100)
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
  #rowData(d_medians)[, "cluster_id"] <-
  #  as.factor(rownames(rowData(d_counts)))

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
  if (is.null(blocking)) {
    res_DA <- testDA_edgeR(d_counts, design, contrast, min_cells = 0)
  } else if (singleCluster) {
    res_DA <- testDA_edgeR(d_counts, design, contrast, min_cells = 0)
  } else {
    daDiagnosticPlotDirectory <- paste0(
      figureDirectory,
      "differentialAbundance",
      clusterName,
      group_id,
      "Visits",
      concatinatedVisits,
      progression,
      siteOfOnset,
      clusterType,
      as.character(singleCluster)
    )

    dir.create(daDiagnosticPlotDirectory, showWarnings = FALSE)

    res_DA <-
      testDA_voom(
        d_counts,
        design,
        contrast,
        block_id = experimentInfo$patient_id,
        min_cells = 0,
        plot = TRUE,
        path = daDiagnosticPlotDirectory
      )
  }

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
      siteOfOnset,
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
      siteOfOnset,
      clusterType,
      "DifferentialAbundanceVolcanoPlot.jpeg"
    )
  )

  # display table of results for top DA clusters
  topTable(res_DA, format_vals = TRUE)

  # calculate number of significant detected DA clusters at 10% false discovery
  # rate (FDR)
  table(topTable(res_DA, all = TRUE)$p_adj <= differentialAbundanceThreshold)

  dsDiagnosticPlotDirectory <- paste0(
    figureDirectory,
    "differentialState",
    clusterName,
    group_id,
    "Visits",
    concatinatedVisits,
    progression,
    siteOfOnset,
    clusterType,
    as.character(singleCluster)
  )

  dir.create(dsDiagnosticPlotDirectory, showWarnings = FALSE)


  # Test for differential states (DS) within clusters
  if (is.null(blocking)) {
    res_DS <-
      testDS_limma(
        d_counts,
        d_medians,
        design,
        contrast,
        min_cells = 0,
        plot = TRUE,
        path = dsDiagnosticPlotDirectory
      )
  } else {
    res_DS <-
      testDS_limma(
        d_counts,
        d_medians,
        design,
        contrast,
        block_id = experimentInfo$patient_id,
        min_cells = 0,
        plot = TRUE,
        path = dsDiagnosticPlotDirectory
      )
  }

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
      siteOfOnset,
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
      siteOfOnset,
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
      siteOfOnset,
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

updateClinicalData <- function(experimentInfo, directoryName){
  if (directoryName == "gpr32BMonocytes" | directoryName == "gpr18Monocytes") {
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
  } else if (directoryName == "gpr32BSenescence" | directoryName == "gpr18Senescence") {
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
  } else if (directoryName == "gpr32TCells") {
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
  } else if (directoryName == "gpr18TCells") {
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS101', arr.ind = TRUE), 'sample_id'] <- 'BAS00101_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS032_02', arr.ind = TRUE), 'sample_id'] <- 'BAS032_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS_033_02', arr.ind = TRUE), 'sample_id'] <- 'BAS033_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS056', arr.ind = TRUE), 'sample_id'] <- 'BAS056_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS_057_02', arr.ind = TRUE), 'sample_id'] <- 'BAS057_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS070', arr.ind = TRUE), 'sample_id'] <- 'BAS070_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS080', arr.ind = TRUE), 'sample_id'] <- 'BAS080_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS083', arr.ind = TRUE), 'sample_id'] <- 'BAS083_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00057_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00057_04_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00057-5', arr.ind = TRUE), 'sample_id'] <- 'BLT00057-5_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00074_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00074_4T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00074-2', arr.ind = TRUE), 'sample_id'] <- 'BLT00074-02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00075-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00075_4T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00092_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00092_04_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00092-6', arr.ind = TRUE), 'sample_id'] <- 'BLT00092-6_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00186-6', arr.ind = TRUE), 'sample_id'] <- 'BLT00186-6_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00186_9', arr.ind = TRUE), 'sample_id'] <- 'BLT00186-9_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00192_02', arr.ind = TRUE), 'sample_id'] <- 'BLT00192_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00198_02', arr.ind = TRUE), 'sample_id'] <- 'BLT00198_02T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00211', arr.ind = TRUE), 'sample_id'] <- 'BLT00211_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00211-3', arr.ind = TRUE), 'sample_id'] <- 'BLT00211-3_ T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00214-2', arr.ind = TRUE), 'sample_id'] <- 'BLT00214-2_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00214_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00214_04'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00214-5', arr.ind = TRUE), 'sample_id'] <- 'BLT00214-5_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00228', arr.ind = TRUE), 'sample_id'] <- 'BLT00228_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00230_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00230_04_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00230', arr.ind = TRUE), 'sample_id'] <- 'BLT00230_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00230-2', arr.ind = TRUE), 'sample_id'] <- 'BLT00230-2_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00243_05', arr.ind = TRUE), 'sample_id'] <- 'BLT00243_05_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00243', arr.ind = TRUE), 'sample_id'] <- 'BLT00243_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00243-7', arr.ind = TRUE), 'sample_id'] <- 'BLT00243-7_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00244_02', arr.ind = TRUE), 'sample_id'] <- 'BLT00244_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00244-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00244-4_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00244-6', arr.ind = TRUE), 'sample_id'] <- 'BLT00244-6_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00254-5', arr.ind = TRUE), 'sample_id'] <- 'BLT00254-5_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00254-7', arr.ind = TRUE), 'sample_id'] <- 'BLT00254-7_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00265', arr.ind = TRUE), 'sample_id'] <- 'BLT00265_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00265-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00265-4_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00265-8', arr.ind = TRUE), 'sample_id'] <- 'BLT00265-8 T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00271-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00271-4_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00271-7', arr.ind = TRUE), 'sample_id'] <- 'BLT00271-7_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00274_02', arr.ind = TRUE), 'sample_id'] <- 'BLT00274_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00274-05', arr.ind = TRUE), 'sample_id'] <- 'BLT00274-5_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00274-6', arr.ind = TRUE), 'sample_id'] <- 'BLT00274-6_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00285-3', arr.ind = TRUE), 'sample_id'] <- 'BLT00285-3_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00285-5', arr.ind = TRUE), 'sample_id'] <- 'BLT00285-5 T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00286', arr.ind = TRUE), 'sample_id'] <- 'BLT00286_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT000286-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00286-4_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00286-6', arr.ind = TRUE), 'sample_id'] <- 'BLT00286-6_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00293', arr.ind = TRUE), 'sample_id'] <- 'BLT00293_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00297_2', arr.ind = TRUE), 'sample_id'] <- 'BLT00297_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00297', arr.ind = TRUE), 'sample_id'] <- 'BLT00297_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00300', arr.ind = TRUE), 'sample_id'] <- 'BLT00300_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00304', arr.ind = TRUE), 'sample_id'] <- 'BLT00304_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00311', arr.ind = TRUE), 'sample_id'] <- 'BLT00311_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BUH00190_03', arr.ind = TRUE), 'sample_id'] <- 'BUH00190_03_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS084', arr.ind = TRUE), 'sample_id'] <- 'HC_BAS084_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00195', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00195_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00209', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00209_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00217', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00217_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00229', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00229_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00236', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00236_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00242_02', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00242_02_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00258', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00258_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00261_2', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00261_2_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00278', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00278_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00284', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00284_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00290', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00290_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00291', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00291_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00296', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00296_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00305', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00305_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00306', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00306_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00307', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00307_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00308', arr.ind = TRUE), 'sample_id'] <- 'HC_BLT00308_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_008', arr.ind = TRUE), 'sample_id'] <- 'QS_008_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_009', arr.ind = TRUE), 'sample_id'] <- 'QS_009_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_014', arr.ind = TRUE), 'sample_id'] <- 'QS_014_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_020', arr.ind = TRUE), 'sample_id'] <- 'QS_020_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_024', arr.ind = TRUE), 'sample_id'] <- 'QS_024_T2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'QS_024-2', arr.ind = TRUE), 'sample_id'] <- 'QS_024-2_T2'
  } else if (directoryName == "gpr18BCells") {
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS_033_02', arr.ind = TRUE), 'sample_id'] <- 'BAS033_02'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BAS_057_02', arr.ind = TRUE), 'sample_id'] <- 'BAS057_02'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00074_04', arr.ind = TRUE), 'sample_id'] <- 'BLT00074-04'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00075-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00075-04'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00186_9', arr.ind = TRUE), 'sample_id'] <- 'BLT00186-9'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00261_2', arr.ind = TRUE), 'sample_id'] <- 'BLT00261_02'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT000286-4', arr.ind = TRUE), 'sample_id'] <- 'BLT00286-04'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00297_2', arr.ind = TRUE), 'sample_id'] <- 'BLT00297-2'
    experimentInfo[which(experimentInfo[, 'sample_id'] == 'BLT00209', arr.ind = TRUE), 'sample_id'] <- 'BLT0209'
  }
  return(experimentInfo)
}

performAllDifferentialAbundanceTests <-
  function(df,
           directoryName,
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

    #'
    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Case vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

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

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Fast Progression vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

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

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = "",
        blocking = NULL
      )
      },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Slow Progression vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

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
        "BulbarLimb",
        "ethnicity",
        "timeFromOnsetToVisitInYears")
    singleCluster <- FALSE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


    ### Fast vs Slow for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Fast Bulbar Onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Fast"
    siteOfOnset <- "Bulbar"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Fast Bulbar Onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Slow Bulbar Onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Slow"
    siteOfOnset <- "Bulbar"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Slow Bulbar Onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Fast Limb Onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Fast"
    siteOfOnset <- "Limb"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Fast Limb Onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Slow Limb Onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- "Slow"
    siteOfOnset <- "Limb"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Slow Limb Onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Bulbar Onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- ""
    siteOfOnset <- "Bulbar"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Bulbar Onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Limb onset vs Controls for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "caseControl"
    visits <- c(1)
    cases <- c("Case", "Control")
    covariants <- c("ageAtVisit", "gender", "ethnicity")
    singleCluster <- FALSE
    progression <- ""
    siteOfOnset <- "Limb"

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Limb onset vs Controls for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = progression, siteOfOnset = siteOfOnset,
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

    ### Bulbar vs Limb for first visit for clusters
    samplesContributionToClustersThreshold <- 10
    differentialAbundanceThreshold <- 0.05
    calculateSampleContributionsToClusters <- FALSE
    group_id <- "BulbarLimb"
    visits <- c(1)
    cases <- c("Case")
    covariants <-
      c("ageAtVisit",
        "gender",
        "fastSlow",
        "ethnicity",
        "timeFromOnsetToVisitInYears")
    singleCluster <- FALSE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


    ### Bulbar vs Limb for first visit for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = NULL
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })#'


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
        "BulbarLimb",
        "timeFromOnsetToVisitInYears",
        "timeFromVisit1InYears"
      )
    singleCluster <- FALSE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
      },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


    ### Visit 1 vs Visit 2 for visit 1 & 2 for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })

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
        "BulbarLimb",
        "timeFromOnsetToVisitInYears",
        "timeFromVisit1InYears"
      )
    singleCluster <- FALSE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
      },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


    ### Visit 2 vs Visit 3 for visit 1 & 3 for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


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
        "BulbarLimb",
        "timeFromOnsetToVisitInYears",
        "timeFromVisit1InYears"
      )
    singleCluster <- FALSE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })


    ### Visit 1 vs Visit 3 for visit 1 & 3 for all cells
    singleCluster <- TRUE

    tryCatch({
      differentialAbundanceAnalysis(
        df = df,
        directoryName = directoryName,
        columnNames = columnNames,
        clusterName = clusterName,
        samplesContributionToClustersThreshold = samplesContributionToClustersThreshold,
        differentialAbundanceThreshold = differentialAbundanceThreshold,
        calculateSampleContributionsToClusters = calculateSampleContributionsToClusters,
        group_id = group_id,
        visits = visits,
        cases = cases,
        covariants = covariants,
        singleCluster = singleCluster,
        markersOrCell = markersOrCell,
        progression = "", siteOfOnset = "",
        blocking = TRUE
      )
    },
    error = function(cond) {
      message(cond)
      setwd("..")
      setwd("..")
      # Choose a return value in case of error
      return(NA)
    })
  }

recalculatePValueAdjustments <-
  function(DA,
           sigCutOff,
           fileNames,
           clusterName,
           markersOrCell,
           directories,
           markerName,
           flipFoldChange = TRUE) {

    names(fileNames) <- c("allCells", "cellPopulations")

    for (directory in directories) {
      i <- 1
      for (file in fileNames) {
        try({
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
            if (directory == "gpr32BCells" | directory == "gpr18BCells") {
              df[, "typeOfCells"] <- "B Cells"
            } else if (directory == "gpr32BMonocytes" | directory == "gpr18Monocytes") {
              df[, "typeOfCells"] <- "Monocytes"
            } else if (directory == "gpr32TCells" | directory == "gpr18TCells") {
              df[, "typeOfCells"] <- "T Cells"
            } else if (directory == "gpr32BSenescence" | directory == "gpr18Senescence") {
              df[, "typeOfCells"] <- "Senescent T Cells"
            }
          } else if (markersOrCell == "Clusters")
          {
            filePath <-
              paste0(
                "data/",
                directory,
                "/clusteringOutput/",
                clusterName,
                "CellPopulations.csv"
              )
            cellPopulations <- read.csv(filePath)

            df <-
              merge(df,
                    cellPopulations[c(clusterName, "cell_population")],
                    by.x = "cluster_id",
                    by.y = clusterName,
                    all.x=TRUE
                    )

            df[, "typeOfCells"] <- df[, "cell_population"]

            df <- df[, -which(names(df) %in% c("cell_population"))]
          } else if (markersOrCell == "Markers")
          {
            filePath <-
              paste0(
                "data/",
                directory,
                "/clusteringOutput/",
                clusterName,
                "CellPopulations.csv"
              )

            cellPopulations <- read.csv(filePath)

            markerFilePath <-
              paste0("data/",
                     directory,
                     "/clusteringOutput/",
                     clusterName,
                     "Markers.csv")

            markerPopulations <- read.csv(markerFilePath)

            colnames(markerPopulations)[ncol(markerPopulations)] <-
              "marker_population"

            cellPopulations <-
              merge(markerPopulations[c(clusterName, "marker_population")],
                    cellPopulations[c(clusterName, "cell_population")],
                    by = clusterName,
                    all.x = TRUE
                    )

            cellPopulations <- cellPopulations[!duplicated(cellPopulations[, c("marker_population", "cell_population")]),]

            df <-
              merge(df,
                    cellPopulations[c("marker_population", "cell_population")],
                    by.x = "cluster_id",
                    by.y = "marker_population",
                    all.x = TRUE
                    )

            df[, "typeOfCells"] <- df[, "cell_population"]

            df <- df[, -which(names(df) %in% c("cell_population"))]
          }
          else {
            df[, "typeOfCells"] <- df[, "cluster_id"]
          }

          variableColumns <- c("AveExpr", "logCPM", "t", "B", "LR")

          emptyColumns <-
            variableColumns[!variableColumns %in% colnames(df)]

          df[, emptyColumns] <- NA

          if (exists("combinedDf")) {
            combinedDf <- rbind(combinedDf, df)
          } else {
            combinedDf <- df
          }
        })
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

    combinedDf$minus_log_p_val <-
      0 - log10(combinedDf[, "p_val"])

    # Define Colours
    mycolors <- data.frame(DOWN = "blue",
                           UP = "red",
                           NO = "black")

    setwd("data")

    dir.create(paste0(markerName, "pValueAdjustmentsResults"), showWarnings = FALSE)
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
      ), quality = 100)
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data = combinedDf,
        aes(
          x = logFC,
          y = minus_log_p_val,
          col = fdr_diff_expressed,
          label = fdr_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential Abundance of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log P-Value")
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
        ), quality = 100
      )
      par(mar = c(1, 1, 1, 1))
      p <- ggplot(
        data =
          combinedDf[combinedDf[, "marker_id"] == "GPR32", ],
        aes(
          x = logFC,
          y = minus_log_p_val,
          col = fdr_diff_expressed,
          label = fdr_label
        )
      ) +
        geom_point() +
        theme_minimal() +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() +
        ggtitle("Differential States of Clusters") +
        xlab("Log Fold Change") + ylab("0 - Log P-Value")
      print(p)
      dev.off()
      gc()
    }

    fwrite(
      combinedDf,
      paste0(
        markerName,
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

calculateClusterMarkers <-
  function(df,
           directoryName,
           clusterName,
           columnNames,
           cutoff,
           markersOrCell = "CellPopulations") {
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

    if (directoryName == "gpr32BCells" | directoryName == "gpr18BCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          IgD_positive == cellPopulationMarkers[cellPopulationMarkersRow, "IgD_positive"] &
            CD24_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD24_positive"] &
            CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32BMonocytes"  | directoryName == "gpr18Monocytes") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          HLA_DR_positive == cellPopulationMarkers[cellPopulationMarkersRow, "HLA_DR_positive"] &
            CD11b_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_positive"] &
            CD11b_activated_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_activated_positive"] &
            CD14_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD14_positive"] &
            CD16_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD16_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32TCells" | directoryName == "gpr18TCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          CD127_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD127_positive"] &
            CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"] &
            CD25_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD25_positive"] &
            FoxP3_positive == cellPopulationMarkers[cellPopulationMarkersRow, "FoxP3_positive"] &
            CD45RO_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RO_positive"] &
            CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32BSenescence" | directoryName == "gpr18Senescence") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        results[results[, clusterName] %in% filter(
          results,
          CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"] &
            CD45RA_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RA_positive"] &
            CD28_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD28_positive"] &
            KLRG1_positive == cellPopulationMarkers[cellPopulationMarkersRow, "KLRG1_positive"] &
            CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"] &
            CCR7_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CCR7_positive"] &
            CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"]
        )[, clusterName], "cell_population"] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    }

    fwrite(
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

    fwrite(
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

    fwrite(
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
      read.csv(paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        "clusteringOutputs.csv"
      ))

    data <- df[, columnNames]
    data[, "ID"] <- row.names(data)
    data <- data[, c("ID", columnNames)]
    head(data)

    clust <- df[, paste0("meta_clusters_flowsom", numberOfClusters)]

    colnames(clust) <- numberOfClusters
    clustColnames <- colnames(clust)
    clust[, "ID"] <- row.names(clust)
    clust[, "1"] <- 1
    clust <- clust[, c("ID", "1", clustColnames)]
    head(clust)

    edist <- function(x) {
      xmean <- mean(x)
      x <- x - xmean
      x <- x ^ 2
      return(x)
    }

    calcVar <- function(clust, data, clust_ID) {
      clust_wIDs <-
        data.frame(clust_ID, clust) #Recombine the clust_ID and clust vectors into a datafram
      colnames(clust_wIDs) <-
        c("ID", "clust") #Rename to appropriate names

      colnames(data)[1] <-
        "ID"               #Rename first column of main dataset to be "ID", since this will be matched within left-join


      #Join the cluster assignments to the main dataset by the ID index
      df <- data %>%
        left_join(clust_wIDs, by = "ID")

      #Compute total variance across all clusters and predictor variables
      totvar <- df %>%
        dplyr::select(-ID) %>%                                           #Drop the ID column
        group_by(clust) %>%                                       #group by clusters
        dplyr::summarise(across(everything(), ~ edist(.))) %>%      #Identify variance for each variable by cluster
        ungroup() %>%                                             #Remove grouping by cluster
        dplyr::select(-clust)                                       #drop cluster identifying column

      print(colSums(totvar))

      totvar <- sum(totvar) #square and sum all variance

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
      ylab("Sum of squared distances") +
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
  function(df,
           clusterName,
           directoryName,
           columnNames,
           markersOrCell,
           markerType,
           prettyColumnNames) {
    # clusterName <- "meta_clusters_flowsom"
    # directoryName <- "gpr32BCells"
    # columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
    # markersOrCell <- "Markers"

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
      merge(cellPopulationMarkers, markerPopulationMarkers, by = clusterName,
            all.x = TRUE)

    populations <-
      populations[, c(clusterName, "cell_population", "marker_population")]

    df <- df[, c(columnNames, clusterName)]

    df[, columnNames] <- scale(df[, columnNames])

    colnames(df)[1:length(columnNames)] <- prettyColumnNames

    colnames(df)[ncol(df)] <- "clusters"

    figureDirectory <- paste0("data/", directoryName, "/figures/")

    df <- df[df[, "clusters"] %in% populations[, clusterName],]
    populations <-
      populations[populations[, clusterName] %in% df[, "clusters"],]

    jpeg(file = paste0(figureDirectory,
                       "na.jpeg"), quality = 100)
    metacluster_colours <-
      as.vector(colortools::wheel('#ce302e', num = nrow(populations)))
    dev.off()

    names(metacluster_colours) <-
      populations[, if (markersOrCell == "CellPopulations") {
        "cell_population"
      } else if (markersOrCell == "Cluster") {
        clusterName
      } else {
        "marker_population"
      }]

    summarisedDf <- df %>%
      group_by(clusters) %>%
      summarise_all(median) %>%
      remove_rownames() %>%
      column_to_rownames('clusters') %>%
      as.matrix() #%>%
    #pheatmap:::scale_rows()

    row_order <-
      populations[match(rownames(summarisedDf), populations[, clusterName]),
                  if (markersOrCell == "CellPopulations") {
                    "cell_population"
                  } else if (markersOrCell == "Cluster") {
                    clusterName
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
        "heatmap",
        clusterName,
        markersOrCell,
        markerType,
        ".jpeg"
      ),
      width = 4000,
      height = 2000,
      res = 320
    )

    print(fig)
    dev.off()
  }

minMaxScaling <- function(x, minValue, maxValue, na.rm = TRUE) {
  return((x - minValue) / (maxValue - minValue))
}

calculateMediansValue <- function(directoryName,
                                  columnNames,
                                  markersOrCell,
                                  clusterName,
                                  df) {
  if (markersOrCell == "CellPopulations") {
    filePath <-
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        "CellPopulations.csv"
      )

    cellPopulationMarkers <- read.csv(filePath)
  } else if (markersOrCell == "Markers") {
    markerFilePath <-
      paste0("data/",
             directoryName,
             "/clusteringOutput/",
             clusterName,
             "Markers.csv")

    cellPopulationMarkers <- read.csv(markerFilePath)
  }
  if (markersOrCell != "Clusters") {
    df <-
      merge(df,
            cellPopulationMarkers[, c(clusterName, "cell_population")],
            by = clusterName,
            all.x = TRUE)
    df[, clusterName] <- df[, "cell_population"]
  }



  columnNamesMedian <- paste0(columnNames, "_median")


  results <-
    data.frame(matrix(ncol = length(columnNamesMedian) + 2,
                      nrow = 0))
  for (filename in unique(df[, "fileName"])) {
    for (cluster in unique(df[, clusterName])) {
      new_row <- c(filename, cluster)

      df2 <-
        df[df[, clusterName] == cluster & df[, "fileName"] == filename, ]

      for (column in columnNames) {
        clusterMedian <- median(df2[, column])

        new_row <- append(new_row, clusterMedian)
      }

      results <- rbind(new_row, results)
    }
  }

  colnames(results) <-
    c("fileName",
      clusterName,
      columnNamesMedian)

  fwrite(
    results,
    paste0(
      "data/",
      directoryName,
      "/clusteringOutput/",
      clusterName,
      markersOrCell,
      "Medians.csv"
    ),
    row.names = FALSE
  )
}

calculateCounts <- function(directoryName,
                                  markersOrCell,
                                  clusterName,
                                  df) {
  if (markersOrCell == "CellPopulations") {
    filePath <-
      paste0(
        "data/",
        directoryName,
        "/clusteringOutput/",
        clusterName,
        "CellPopulations.csv"
      )

    cellPopulationMarkers <- read.csv(filePath)
  } else if (markersOrCell == "Markers") {
    markerFilePath <-
      paste0("data/",
             directoryName,
             "/clusteringOutput/",
             clusterName,
             "Markers.csv")

    cellPopulationMarkers <- read.csv(markerFilePath)
  }
  if (markersOrCell != "Clusters") {
    df <-
      merge(df,
            cellPopulationMarkers[, c(clusterName, "cell_population")],
            by = clusterName,
            all.x = TRUE)
    df[, clusterName] <- df[, "cell_population"]
  }

  results <-
    data.frame(matrix(ncol = 3,
                      nrow = 0))

  for (filename in unique(df[, "fileName"])) {
    for (cluster in unique(df[, clusterName])) {
      df2 <-
        df[df[, clusterName] == cluster & df[, "fileName"] == filename, ]

      count <- nrow(df2)
      new_row <- c(filename, cluster, count)

      results <- rbind(results, new_row)
    }
  }

  colnames(results) <-
    c("fileName",
      clusterName,
      "count")

  fwrite(
    results,
    paste0(
      "data/",
      directoryName,
      "/clusteringOutput/",
      clusterName,
      markersOrCell,
      "Counts.csv"
    ),
    row.names = FALSE
  )
}

identifyUnstableClustersFromCounts <- function(directoryName,
                            markersOrCell,
                            clusterName,
                            df, cutOff) {
  results <-
    data.frame(matrix(ncol = 3,
                      nrow = 0))

  colnames(results) <-
    c("fileName",
      clusterName,
      "count")

  df[, "count"] <- as.double(df[, "count"])

  for (cluster in unique(df[, clusterName])) {
    df2 <- df[df[, clusterName] ==  cluster,]
    totalCount <- sum(df2[, "count"])

    for (filename in unique(df2[, "fileName"])) {
      df3 <- df2[df2[, "fileName"] == filename, ]
      if (df3[, "count"] > cutOff * totalCount) {
        results <- rbind(results, df3)
      }
    }
  }

  fwrite(
    results,
    paste0(
      "data/",
      directoryName,
      "/clusteringOutput/",
      clusterName,
      markersOrCell,
      "CountsOverCutoff.csv"
    ),
    row.names = FALSE
  )
}


differentialCombinedManhattanPlot <- function(pattern, clusterName, figureName, markersOrCell, directoryNames, markerName) {
  directory <- paste0("data/", markerName, "pValueAdjustmentsResults/")

  fileNames <- list.files(directory, pattern = pattern)

  fileNames <- fileNames[grepl(clusterName, fileNames, fixed = TRUE)]

  fileNamesPath <- paste0(directory, fileNames)

  dfs <- lapply(fileNamesPath, read.csv)

  names(dfs) <- fileNames

  i <- 1

  for(df in dfs){
    df$experiment <- names(dfs)[i]
    i <- i + 1

    if (exists("combinedDf")){
      combinedDf <- rbind(combinedDf, df)
    } else {
      combinedDf <- df
    }
  }

  for(directoryName in directoryNames){
    stability <- fread(paste0("data/", directoryName, "/clusteringOutput/", clusterName, "MarkerStability", ".csv"))
    stability <- as.data.frame(stability)
    numberOfColumns <- ncol(stability)

    stability$median <- NA

    for(row in seq(nrow(stability))) {
      stability[row, "median"] <- median(as.double(stability[row, seq(from=numberOfColumns-99, to=numberOfColumns)]))
    }

    clustermarkerName <- paste0(clusterName, "Marker")

    if(!exists("combinedStability")){
      combinedStability <- stability[,c(clustermarkerName, "median")]
    } else {
      combinedStability <- rbind(combinedStability, stability[,c(clustermarkerName, "median")])
    }
  }

  combinedDf <- merge(combinedDf, combinedStability, by.x = "cluster_id", by.y = clustermarkerName, all.x = TRUE)

  combinedDf[combinedDf$cluster_id == 1, "median"] <- 1

  combinedDf$jaccard <- NA
  combinedDf$shape <- NA

  combinedDf[is.na(combinedDf$median), "median"] <- 0

  combinedDf[combinedDf$median < 0.65, "jaccard"] <- "Low"
  combinedDf[combinedDf$median > 0.65, "jaccard"] <- "Moderate"
  combinedDf[combinedDf$median > 0.75, "jaccard"] <- "High"
  combinedDf[combinedDf$median > 0.85, "jaccard"] <- "Very High"

  combinedDf[combinedDf$jaccard == "Low", "shape"] <- 18
  combinedDf[combinedDf$jaccard == "Moderate", "shape"] <- 16
  combinedDf[combinedDf$jaccard == "High", "shape"] <- 15
  combinedDf[combinedDf$jaccard == "Very High", "shape"] <- 17

  combinedDf <- combinedDf[order(combinedDf$median),]


  combinedDf$jaccard <- factor(combinedDf$jaccard, levels = c("Very High", "High", "Moderate", "Low"))

  combinedDf <- combinedDf[combinedDf$typeOfCells %in% combinedDf[combinedDf$fdr_adjusted_p_val<0.05, "typeOfCells"], ]
  combinedDf <- combinedDf[combinedDf$cluster_id %in% combinedDf[combinedDf$fdr_adjusted_p_val<0.05, "cluster_id"], ]
  combinedDf <- combinedDf[combinedDf$experiment %in% combinedDf[combinedDf$fdr_adjusted_p_val<0.05, "experiment"], ]
  # combinedDf[combinedDf$cluster_id == "CD8+ CD4- CD27+ CD45RA+ KLRG1- CCR7+ CD28- T Cells", "typeOfCells"] <- "Naïve CD8+ T Cells (Senescence Panel)"
  # combinedDf[combinedDf$cluster_id == "CD8+ CD4- CD45RO- CD25- CD127- FoxP3+ CD3+ T Cells", "typeOfCells"] <- "Naïve CD8+ T Cells"
  # combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25+ CD127+ FoxP3+ CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (FoxP3+)"
  # combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25+ CD127+ FoxP3- CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (FoxP3-)"
  # combinedDf[combinedDf$cluster_id == "CD27+ CD24- IgD- CD19+ B Cells", "typeOfCells"] <- "Switched Memory B Cells (CD24-)"
  # combinedDf[combinedDf$cluster_id == "CD27+ CD24+ IgD- CD19+ B Cells", "typeOfCells"] <- "Switched Memory B Cells (CD24+)"
  # combinedDf[combinedDf$cluster_id == "CD27+ CD24- IgD+ CD19+ B Cells", "typeOfCells"] <- "Unswitched Memory B Cells (CD24-)"
  # combinedDf[combinedDf$cluster_id == "CD27+ CD24+ IgD+ CD19+ B Cells", "typeOfCells"] <- "Unswitched Memory B Cells (CD24+)"
  # combinedDf[combinedDf$cluster_id == "HLA-DR-/low CD16+ CD14+ CD11b- CD11b Activated+ Monocytes", "typeOfCells"] <- "HLA-DR Negative/Low Activated CD11b+ Intermediate Monocytes (CD11b Low)"
  # combinedDf[combinedDf$cluster_id == "HLA-DR-/low CD16+ CD14+ CD11b+ CD11b Activated+ Monocytes", "typeOfCells"] <- "HLA-DR Negative/Low Activated CD11b+ Intermediate Monocytes (CD11b High)"
  # combinedDf[combinedDf$cluster_id == "CD8+ CD4- CD27+ CD45RA+ KLRG1- CCR7+ CD28+ T Cells", "typeOfCells"] <- "Naïve CD8+ T Cells (Senescence Panel)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25+ CD127+ FoxP3- CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25+ CD127+ FoxP3-)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25+ CD127+ FoxP3+ CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25+ CD127+ FoxP3+)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25- CD127+ FoxP3+ CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25- CD127+ FoxP3+)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25- CD127- FoxP3+ CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25- CD127- FoxP3+)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25- CD127- FoxP3- CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25- CD127- FoxP3-)"
  combinedDf[combinedDf$cluster_id == "CD8- CD4- CD45RO+ CD25- CD127- FoxP3- CD3+ T Cells", "typeOfCells"] <- "Memory Double Negative T Cells (CD25- CD127- FoxP3-)"

  # message(unique(combinedDf$cluster_id))

  combinedDf$typeOfCells <- gsub(" CD11b\\+", " CD11b+\n", combinedDf$typeOfCells)
  combinedDf$typeOfCells <- gsub("Memory Double Negative T Cells ", "Memory Double Negative T Cells\n", combinedDf$typeOfCells)
  combinedDf$typeOfCells <- gsub("B Cells ", "B Cells\n", combinedDf$typeOfCells)
  combinedDf$typeOfCells <- gsub("Memory Double Negative Regulatory ", "Memory Double Negative Regulatory\n", combinedDf$typeOfCells)

  unique(combinedDf[combinedDf$cluster_id %in% c("CD8+ CD4- CD27+ CD45RA+ KLRG1- CCR7+ CD28+ T Cells",
                               "CD8- CD4+ CD45RO- CD25- CD127- FoxP3- CD3+ T Cells",
                               "CD8- CD4+ CD45RO+ CD25- CD127+ FoxP3+ CD3+ T Cells",
                               "CD8- CD4- CD45RO+ CD25- CD127- FoxP3+ CD3+ T Cells",
                               "CD8+ CD4- CD45RO+ CD25+ CD127+ FoxP3- CD3+ T Cells",
                               "CD8- CD4- CD45RO+ CD25+ CD127+ FoxP3- CD3+ T Cells",
                               "CD8+ CD4- CD45RO+ CD25+ CD127+ FoxP3+ CD3+ T Cells",
                               "CD8- CD4+ CD45RO+ CD25+ CD127- FoxP3+ CD3+ T Cells",
                               "CD8- CD4- CD45RO- CD25- CD127- FoxP3+ CD3+ T Cells"),"experiment"])

  #unique(combinedDf[order(combinedDf$typeOfCells), c("cluster_id", "panel", "typeOfCells")])

  combinedDf[combinedDf$experiment == paste0(clusterName,"BulbarLimbVisits1AllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Bulbar vs Limb"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1AllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "ALS vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1BulbarAllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Bulbar vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastAllCells", figureName, markersOrCell, ".csv"), "experiment"] <- "Fast vs Controls"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastBulbarAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast Bulbar vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1FastLimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast Limb vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1LimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Limb vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowBulbarAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow Bulbar vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"caseControlVisits1SlowLimbAllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Slow Limb vs Control"
  combinedDf[combinedDf$experiment == paste0(clusterName,"fastSlowVisits1AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Fast vs Slow"
  combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits12AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 2 vs Visit 1"
  combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits13AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 3 vs Visit 1"
  combinedDf[combinedDf$experiment == paste0(clusterName,"visitVisits23AllCells", figureName, markersOrCell, ".csv"),"experiment"] <- "Visit 2 vs Visit 3"

  combinedDf$panel <- factor(combinedDf$panel)
  combinedDf <- combinedDf[order(combinedDf$typeOfCells),]
  combinedDf <- combinedDf[order(combinedDf$panel),]

  combinedDf$typeOfCells <- factor(combinedDf$typeOfCells, levels =  unique(combinedDf$typeOfCells))

  #combinedDf[combinedDf$typeOfCells == "Memory Double Negative T Cells", ]

  #combinedDf$typeOfCells <- unlist(lapply(combinedDf$typeOfCells, stringBreak, sep = " ", buffer = 50))

  #dir.create("data/combinedFigures", showWarnings = FALSE)
  #jpeg(filename = paste0("data/combinedFigures/", clusterName, figureName, markersOrCell, ".jpeg"))
  fig <-
    ggplot(
      combinedDf,
      aes(
        x = as.factor(typeOfCells),
        y = minus_log_fdr_adjusted_p_val,
        color = logFC,
        shape = jaccard
      )
    ) +
    geom_point(alpha = 0.75, size=5) +
    scale_shape_manual(values=unique(combinedDf$shape)) +
    facet_wrap(~experiment) +
    guides(color = guide_colourbar(title="log2(Fold Change)", order = 1)) +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 10
    )) +
    xlab("Cell Populations") +
    ylab("-log10(Adjusted P-Value)") +
    ylim(0, NA) + guides(
      shape = guide_legend(title = "Stability")
    ) +
    geom_hline(yintercept = 0 - log10(0.05), linetype = "dashed") +
    scale_colour_viridis_c()

  #dev.off()
  return(fig)
}


generateSubsampledPhenographClusters <- function(directoryName,
                                                 columnNames,
                                                 clusterName,
                                                 knn) {
  tryCatch({
    d_f <<-
      fread(file=
              paste0(
                "data/",
                directoryName,
                "/clusteringOutput/phenographClusterStability.csv"
              )
      )
    d_f <<- as.data.frame(d_f)

    iterationsMin <<-
      max(na.omit(suppressWarnings(as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(d_f)[!colnames(d_f) %in% columnNames]))))) + 1
    iterationsMax <<- iterationsMin

  }, warning = function(cond) {
    d_f <<-
      fread(file=paste0(
        "data/",
        directoryName,
        "/clusteringOutput/clusteringOutputs.csv"
      ))
    d_f <<- as.data.frame(d_f)

    d_f <<- d_f[, c("fileName", columnNames, clusterName)]

    iterationsMin <<- 1
    iterationsMax <<- 1

  }, error = function(cond) {
    d_f <<-
      fread(file=paste0(
        "data/",
        directoryName,
        "/clusteringOutput/clusteringOutputs.csv"
      ))
    d_f <<- as.data.frame(d_f)

    d_f <<- d_f[, c("fileName", columnNames, clusterName)]

    iterationsMin <<- 1
    iterationsMax <<- 1
  })

  message(paste0("Max iterations:", iterationsMax))

  experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

  experimentInfo <- as.data.frame(experimentInfo)

  experimentInfo <- experimentInfo[experimentInfo$experiment == "flowCytometry", ]

  experimentInfo <- experimentInfo[, c("patient_id", "sample_id")]

  experimentInfo <- updateClinicalData(experimentInfo, directoryName)

  fileNames <- unique(d_f[, "fileName"])

  experimentInfo <- experimentInfo[experimentInfo$sample_id %in% fileNames,]

  # my.cluster <- parallel::makeCluster(n.cores)
  # doParallel::registerDoParallel(cl = my.cluster)
  # message("Server cluster defined")
  # message(foreach::getDoParRegistered())
  # message(foreach::getDoParWorkers())
  #
  # clusterResults <-
  #  foreach(number = seq(from = iterationsMin, to = iterationsMax),
  #          .combine = 'cbind',
  #          .errorhandling = "remove") %dopar% {
  for (number in seq(from = iterationsMin, to = iterationsMax)) {
    try({
      # try(source("R/01_functions.R"))
      # try(source("R/00_datasets.R"))
      # loadlibraries()

      randomNumbers <-
        sample(seq(length(unique(
          experimentInfo$patient_id
        ))), floor(0.2 * length(unique(
          experimentInfo$patient_id
        ))), replace =
          FALSE)

      message("Random Numbers")
      message(randomNumbers)
      length(randomNumbers)

      keptFileNames <-
        experimentInfo[!experimentInfo$sample_id %in% experimentInfo[experimentInfo$patient_id %in% unique(experimentInfo[, "patient_id"])[randomNumbers], "sample_id"], "sample_id"]

      df2 <- d_f[d_f[, "fileName"] %in% keptFileNames,]

      gc(full = TRUE)

      phenograph <-
        Rphenograph(df2[, columnNames], k = knn)

      ## attributes(phenograph)

      clusters_phenograph <- phenograph$membership

      if(iterationsMin > 1) {
        d_f <<-
          fread(file=
                  paste0(
                    "data/",
                    directoryName,
                    "/clusteringOutput/phenographClusterStability.csv"
                  )
          )

        d_f <<- as.data.frame(d_f)

        clusterNumber <-
          max(na.omit(suppressWarnings(as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(d_f)[!colnames(d_f) %in% columnNames]))))) + 1

      } else {
        tryCatch({
          d_f <<-
            fread(file=
                    paste0(
                      "data/",
                      directoryName,
                      "/clusteringOutput/phenographClusterStability.csv"
                    )
            )

          d_f <<- as.data.frame(d_f)

          clusterNumber <<-
            max(na.omit(suppressWarnings(as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(d_f)[!colnames(d_f) %in% columnNames]))))) + 1

        },
        error=function(cond) {
          clusterNumber <<- number
        })
      }

      d_f[, paste0(clusterName, clusterNumber)] <- NA

      try({
        d_f[d_f[, "fileName"] %in% keptFileNames, paste0(clusterName, clusterNumber)] <-
          clusters_phenograph
      })

      fwrite(
        d_f,
        paste0(
          "data/",
          directoryName,
          "/clusteringOutput/phenographClusterStability.csv"
        ),
        row.names = FALSE
      )

      #return(d_f[, paste0(clusterName, number)])
    })
  }

  # clusterResults <- as.data.frame(clusterResults)
  # try(colnames(clusterResults) <- paste0(clusterName, seq(from = iterationsMin, to = ncol(clusterResults) - 1 + iterationsMin)))
  # d_f <- cbind(d_f, clusterResults)

}


generateSubsampledFlowsomClusters <- function(directoryName,
                                               columnNames,
                                               clusterNames,
                                               numberOfClusters,
                                               iterations) {

  experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

  experimentInfo <- as.data.frame(experimentInfo)

  experimentInfo <- experimentInfo[experimentInfo$experiment == "flowCytometry", ]

  experimentInfo <- experimentInfo[, c("patient_id", "sample_id")]

  experimentInfo <- updateClinicalData(experimentInfo, directoryName)

  experimentInfo$sample_id <- paste0(experimentInfo$sample_id, ".fcs")

  df <-
    fread(paste0(
      "data/",
      directoryName,
      "/clusteringOutput/clusteringOutputs.csv"
    ))

  df <- as.data.frame(df)

  df <- df[, c("fileName", columnNames, clusterNames)]

  workingDirectory <- getwd()

  setwd(paste0("./data/", directoryName))

  dirFCS <- paste0(getwd(), "/dataPPOutput/scaledFcs")

  tempDirectory <- paste0(getwd(), "/dataPPOutput/temp/")

  fileNames <- list.files(dirFCS)

  fileNames <- fileNames[fileNames != "annotation.txt"]

  experimentInfo <- experimentInfo[experimentInfo$sample_id %in% fileNames,]

  message("Test If FileNames and Experiment FileNames Match")
  message(length(experimentInfo$sample_id) == length(fileNames))
  message(all(fileNames %in% experimentInfo$sample_id))
  message(all(experimentInfo$sample_id %in% fileNames))

  dir.create(tempDirectory, showWarnings = FALSE)

  columnIndexes <- seq(length(columnNames))

  # for (number in seq(from=91, to = iterations))

  for (number in seq(iterations)) {
    message(paste0("Iteration :", number))
    try({
      try(file.remove(paste0(tempDirectory, "/", experimentInfo$sample_id)))

      file.copy(
        from = paste0(dirFCS, "/", experimentInfo$sample_id),
        to = tempDirectory,
        overwrite = TRUE,
        recursive = TRUE,
        copy.mode = TRUE
      )

      randomNumbers <-
        sample(seq(length(unique(
          experimentInfo$patient_id
        ))), floor(0.2 * length(unique(
          experimentInfo$patient_id
        ))), replace =
          FALSE)

      message(randomNumbers)

      removedFileNames <-
        experimentInfo[experimentInfo$sample_id %in% experimentInfo[experimentInfo$patient_id %in% unique(experimentInfo[, "patient_id"])[randomNumbers], "sample_id"], "sample_id"]
      keptFileNames <-
        experimentInfo[!experimentInfo$sample_id %in% experimentInfo[experimentInfo$patient_id %in% unique(experimentInfo[, "patient_id"])[randomNumbers], "sample_id"], "sample_id"]

      message("Check No Overlap between Kept and Removed Files")
      message(!all(keptFileNames %in% removedFileNames))
      message(!all(removedFileNames %in% keptFileNames))

      try(file.remove(paste0(tempDirectory, "/", removedFileNames)))

      keptFileNames <- str_replace_all(keptFileNames, ".fcs", "")

      #message(keptFileNames)

      df[, paste0(clusterNames, number)] <- NA

      #run flowsom
      flowsom <- FlowSOM(
        input = tempDirectory,
        transform = FALSE,
        scale = FALSE,
        colsToUse = columnIndexes,
        #seed = seed,
        nClus = numberOfClusters
      )

      # Get metaclustering per cell
      clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
      meta_clusters_flowsom <- as.factor(flowsom$map$mapping[, 1])
      levels(meta_clusters_flowsom) <- flowsom$metaclustering

      df[df[, "fileName"] %in% keptFileNames, paste0(clusterNames, number)] <-
        c(as.integer(clusters_flowsom),
          as.integer(meta_clusters_flowsom))

      if (number %% 10 == 0) {
        fwrite(df,
                  'clusteringOutput/flowsomClusterStability.csv',
                  row.names = FALSE)
      }
    })
  }

  setwd(workingDirectory)
}

identifyPhenographClusterSimilarity <- function(df,
                                             results,
                                             directoryName,
                                             clusterName,
                                             columnNames) {
  df <- df[, !colSums(is.na(df)) == nrow(df)]

  if (clusterName == "clusters_phenographMarker") {
    cellPopulationMarkers <-
      fread(file=
        paste0(
          "./data/",
          directoryName,
          '/clusteringOutput/',
          "clusters_phenograph",
          "Markers",
          '.csv'
        )
      )

    cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)

    df <-
      merge(df, cellPopulationMarkers[, c("clusters_phenograph", "cell_population")], by = "clusters_phenograph")

    colnames(df)[colnames(df) == "cell_population"] <-
      clusterName

  }

  iterations <-
    na.omit(suppressWarnings(as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(df)[!colnames(df) %in% columnNames]))))

  columnOfnterestName <- clusterName

  resultsDf <-
    data.frame(row.names = unique(df[, columnOfnterestName]))

  resultsDf[, columnOfnterestName] <-
    unique(df[, columnOfnterestName])


  for (number in iterations) {
    try({
      message(paste0("\nIteration: ", number))
      column <- paste0(columnOfnterestName, number)

      resultsDf[, column] <- NA

      if (columnOfnterestName != "clusters_phenograph") {
        filteredDf <-
          merge(
            df,
            results[, c("clusters_phenograph", paste0(c("clusters_phenographMarker"),
                                                      number))],
            by.x = paste0("clusters_phenograph", number),
            by.y = "clusters_phenograph"
          )
      } else {
        filteredDf <- df
      }

      filteredDf <-
        na.omit(filteredDf[, c(columnOfnterestName, column)])

      tabulatedDf <- table(filteredDf)
      tabulatedDf <- as.data.frame(tabulatedDf)

      for (cluster in unique(df[, columnOfnterestName])) {
        message(paste0("Cluster: ", cluster))
        filteredTabulatedDf <-
          tabulatedDf[tabulatedDf[, columnOfnterestName] == cluster,]

        jaccardIndexList <- c()

        for (boostaqppedCluster in unique(filteredTabulatedDf[, column])) {
          numerator <-
            filteredTabulatedDf[filteredTabulatedDf[, column] == boostaqppedCluster, "Freq"]
          denominator <-
            sum(filteredTabulatedDf[, "Freq"]) + sum(tabulatedDf[tabulatedDf[, column] == boostaqppedCluster, "Freq"]) - numerator

          jaccardIndex <- numerator / denominator

          jaccardIndexList <- append(jaccardIndexList, jaccardIndex)
        }

        resultsDf[resultsDf[, columnOfnterestName] == cluster, column] <-
          max(jaccardIndexList)
      }
    })
  }

  fwrite(
    resultsDf,
    paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/',
      clusterName,
      "Stability.csv"
    ),
    row.names = FALSE
  )
}

identifyFlowsomClusterSimilarity <- function(df,
                                             results,
                                             directoryName,
                                             clusterName,
                                             iterations) {
  if (clusterName == "meta_clusters_flowsomMarker") {
    cellPopulationMarkers <-
      fread(file=
        paste0(
          "./data/",
          directoryName,
          '/clusteringOutput/',
          "meta_clusters_flowsom",
          "Markers",
          '.csv'
        )
      )

    cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)

    df <-
      merge(df, cellPopulationMarkers[, c("meta_clusters_flowsom", "cell_population")], by = "meta_clusters_flowsom")

    colnames(df)[colnames(df) == "cell_population"] <-
      clusterName

  }

  columnOfnterestName <- clusterName

  resultsDf <-
    data.frame(row.names = unique(df[, columnOfnterestName]))

  resultsDf[, columnOfnterestName] <-
    unique(df[, columnOfnterestName])

  try({
    for (number in seq(iterations)) {
      message(paste0("\nIteration: ", number))
      column <- paste0(columnOfnterestName, number)

      resultsDf[, column] <- NA

      if (columnOfnterestName != "clusters_flowsom") {
        filteredDf <-
          merge(df,
                results[, c("clusters_flowsom", paste0(
                  c(
                    "meta_clusters_flowsom",
                    "meta_clusters_flowsomMarker"
                  ),
                  number
                ))],
                by.x = paste0("clusters_flowsom", number),
                by.y = "clusters_flowsom")
      } else {
        filteredDf <- df
      }

      filteredDf <-
        na.omit(filteredDf[, c(columnOfnterestName, column)])

      tabulatedDf <- table(filteredDf)
      tabulatedDf <- as.data.frame(tabulatedDf)

      for (cluster in unique(df[, columnOfnterestName])) {
        message(paste0("Cluster: ", cluster))
        filteredTabulatedDf <-
          tabulatedDf[tabulatedDf[, columnOfnterestName] == cluster, ]

        jaccardIndexList <- c()

        for (boostaqppedCluster in unique(filteredTabulatedDf[, column])){
          numerator <- filteredTabulatedDf[filteredTabulatedDf[, column] == boostaqppedCluster, "Freq"]
          denominator <- sum(filteredTabulatedDf[, "Freq"]) + sum(tabulatedDf[tabulatedDf[, column] == boostaqppedCluster, "Freq"]) - numerator

          jaccardIndex <- numerator/denominator

          jaccardIndexList <- append(jaccardIndexList, jaccardIndex)
        }

        resultsDf[resultsDf[, columnOfnterestName] == cluster, column] <-
          max(jaccardIndexList)
      }
    }
  })

  fwrite(
    resultsDf,
    paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/',
      clusterName,
      "Stability.csv"
    ),
    row.names = FALSE
  )


}

identifyPhenographBoostrappedCellPopulations <- function(df,
                                                      directoryName,
                                                      iterations,
                                                      markersOrCell,
                                                      cutoff,
                                                      columnNames){
  df <- df[, !colSums(is.na(df)) == nrow(df)]

  columnNamesMedian <- paste0(columnNames, "_median")
  columnNamesPositive <- paste0(columnNames, "_positive")

  iterations <<-
    na.omit(suppressWarnings(as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(df)[!colnames(df) %in% columnNames]))))

  for (number in iterations) {
    try({
      message("")
      message(paste0("Iteration:", number))
      clusterName <- paste0("clusters_phenograph", number)
      markerName <- paste0("clusters_phenographMarker", number)

      subsampledDf <-
        na.omit(df[, c(columnNames, "clusters_phenograph", clusterName)])

      markerResults <-
        data.frame(matrix(
          ncol = 2 * length(columnNamesMedian) + 2,
          nrow = 0
        ))

      for (cluster in unique(subsampledDf[, clusterName])) {
        message(paste0("Cluster:", cluster))

        new_row <- c(cluster)

        df2 <-
          subsampledDf[subsampledDf[, clusterName] == cluster, ]

        for (column in columnNames) {
          clusterMedian <- median(na.omit(df2[, column]))

          new_row <- append(new_row, clusterMedian)
        }

        new_row <-
          append(new_row, replicate(length(colnames(
            markerResults
          )) - length(new_row), NA))

        markerResults <- rbind(new_row, markerResults)
      }

      colnames(markerResults) <-
        c(clusterName,
          columnNamesMedian,
          columnNamesPositive,
          markersOrCell)

      i <- 1

      for (column in columnNamesMedian) {
        markerResults[, columnNamesPositive[i]] <-
          markerResults[, columnNamesMedian[i]] > cutoff[i]
        i <- i + 1
      }

      if (markersOrCell == "CellPopulations") {
        cellPopulationMarkers <-
          fread(file = paste0("data/metadata/", directoryName, ".csv"))
        cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)
      } else if (markersOrCell == "Markers") {
        cellPopulationMarkers <-
          fread(file = paste0("data/metadata/", directoryName, "Markers.csv"))
        cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)
      }

      cellPopulationMarkers <-
        cellPopulationMarkers[, c("name", columnNamesPositive)]

      if (directoryName == "gpr32BCells" | directoryName == "gpr18BCells") {
        for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
          markerResults[markerResults[, clusterName] %in% filter(
            markerResults,
            IgD_positive == cellPopulationMarkers[cellPopulationMarkersRow, "IgD_positive"] &
              CD24_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD24_positive"] &
              CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"]
          )[, clusterName], markerName] <-
            cellPopulationMarkers[cellPopulationMarkersRow, "name"]
        }
      } else if (directoryName == "gpr32BMonocytes" | directoryName == "gpr18Monocytes") {
        for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
          markerResults[markerResults[, clusterName] %in% filter(
            markerResults,
            HLA_DR_positive == cellPopulationMarkers[cellPopulationMarkersRow, "HLA_DR_positive"] &
              CD11b_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_positive"] &
              CD11b_activated_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_activated_positive"] &
              CD14_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD14_positive"] &
              CD16_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD16_positive"]
          )[, clusterName], markerName] <-
            cellPopulationMarkers[cellPopulationMarkersRow, "name"]
        }
      } else if (directoryName == "gpr32TCells" | directoryName == "gpr18TCells") {
        for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
          markerResults[markerResults[, clusterName] %in% filter(
            markerResults,
            CD127_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD127_positive"] &
              CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"] &
              CD25_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD25_positive"] &
              FoxP3_positive == cellPopulationMarkers[cellPopulationMarkersRow, "FoxP3_positive"] &
              CD45RO_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RO_positive"] &
              CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"]
          )[, clusterName], markerName] <-
            cellPopulationMarkers[cellPopulationMarkersRow, "name"]
        }
      } else if (directoryName == "gpr32BSenescence" | directoryName == "gpr18Senescence") {
        for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
          markerResults[markerResults[, clusterName] %in% filter(
            markerResults,
            CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"] &
              CD45RA_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RA_positive"] &
              CD28_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD28_positive"] &
              KLRG1_positive == cellPopulationMarkers[cellPopulationMarkersRow, "KLRG1_positive"] &
              CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"] &
              CCR7_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CCR7_positive"] &
              CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"]
          )[, clusterName], markerName] <-
            cellPopulationMarkers[cellPopulationMarkersRow, "name"]
        }
      }

      if (exists("results")) {
        results <-
          merge(
            results,
            markerResults[, c(clusterName, markerName)],
            by.x = "clusters_phenograph",
            by.y = clusterName,
            all.x = TRUE
          )
      } else {
        results <- seq(1:100)
        results <- as.data.frame(results)
        colnames(results) <- "clusters_phenograph"

        results <-
          merge(
            results,
            markerResults[, c(clusterName, markerName)],
            by.x = "clusters_phenograph",
            by.y = clusterName,
            all.x = TRUE
          )
      }
    })
  }
  fwrite(
    results,
    paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/',
      'clusters_phenograph_Stability.csv'
    ),
    row.names =  FALSE
  )
}

identifyFlowSomBoostrappedCellPopulations <- function(df,
                                                      directoryName,
                                                      iterations,
                                                      markersOrCell,
                                                      cutoff,
                                                      columnNames){
  results <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/', 'meta_clusters_flowsom_Stability.csv'))
  results <- as.data.frame(results)

  columnNamesMedian <- paste0(columnNames, "_median")
  columnNamesPositive <- paste0(columnNames, "_positive")

  for (number in seq(iterations)) {
    message("")
    message(paste0("Iteration:", number))
    metaClustername <- paste0("meta_clusters_flowsom", number)
    clusterName <- paste0("clusters_flowsom", number)
    markerName <- paste0("meta_clusters_flowsomMarker", number)

    subsampledDf <-
      na.omit(df[, c(columnNames, "clusters_flowsom", clusterName)])

    subsampledDf <-
      merge(subsampledDf, results[, c("clusters_flowsom", metaClustername)], by.x = clusterName, by.y = "clusters_flowsom")

    markerResults <-
      data.frame(matrix(
        ncol = 2 * length(columnNamesMedian) + 2,
        nrow = 0
      ))

    for (cluster in unique(subsampledDf[, metaClustername])) {
      message(paste0("Cluster:", cluster))

      new_row <- c(cluster)

      df2 <-
        subsampledDf[subsampledDf[, metaClustername] == cluster,]

      for (column in columnNames) {
        clusterMedian <- median(na.omit(df2[, column]))

        new_row <- append(new_row, clusterMedian)
      }

      new_row <-
        append(new_row, replicate(length(colnames(markerResults)) - length(new_row), NA))

      markerResults <- rbind(new_row, markerResults)
    }

    colnames(markerResults) <-
      c(metaClustername,
        columnNamesMedian,
        columnNamesPositive,
        markersOrCell)

    i <- 1

    for (column in columnNamesMedian) {
      markerResults[, columnNamesPositive[i]] <-
        markerResults[, columnNamesMedian[i]] > cutoff[i]
      i <- i + 1
    }

    if (markersOrCell == "CellPopulations") {
      cellPopulationMarkers <-
        fread(file=paste0("data/metadata/", directoryName, ".csv"))
      cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)
    } else if (markersOrCell == "Markers") {
      cellPopulationMarkers <-
        fread(file=paste0("data/metadata/", directoryName, "Markers.csv"))
      cellPopulationMarkers <- as.data.frame(cellPopulationMarkers)
      }

    cellPopulationMarkers <-
      cellPopulationMarkers[, c("name", columnNamesPositive)]

    if (directoryName == "gpr32BCells" | directoryName == "gpr18BCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        markerResults[markerResults[, metaClustername] %in% filter(
          markerResults,
          IgD_positive == cellPopulationMarkers[cellPopulationMarkersRow, "IgD_positive"] &
            CD24_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD24_positive"] &
            CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"]
        )[, metaClustername], markerName] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32BMonocytes" | directoryName == "gpr18Monocytes") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        markerResults[markerResults[, metaClustername] %in% filter(
          markerResults,
          HLA_DR_positive == cellPopulationMarkers[cellPopulationMarkersRow, "HLA_DR_positive"] &
            CD11b_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_positive"] &
            CD11b_activated_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD11b_activated_positive"] &
            CD14_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD14_positive"] &
            CD16_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD16_positive"]
        )[, metaClustername], markerName] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32TCells" | directoryName == "gpr18TCells") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        markerResults[markerResults[, metaClustername] %in% filter(
          markerResults,
          CD127_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD127_positive"] &
            CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"] &
            CD25_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD25_positive"] &
            FoxP3_positive == cellPopulationMarkers[cellPopulationMarkersRow, "FoxP3_positive"] &
            CD45RO_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RO_positive"] &
            CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"]
        )[, metaClustername], markerName] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    } else if (directoryName == "gpr32BSenescence" | directoryName == "gpr18Senescence") {
      for (cellPopulationMarkersRow in seq(nrow(cellPopulationMarkers))) {
        markerResults[markerResults[, metaClustername] %in% filter(
          markerResults,
          CD27_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD27_positive"] &
            CD45RA_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD45RA_positive"] &
            CD28_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD28_positive"] &
            KLRG1_positive == cellPopulationMarkers[cellPopulationMarkersRow, "KLRG1_positive"] &
            CD4_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD4_positive"] &
            CCR7_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CCR7_positive"] &
            CD8_positive == cellPopulationMarkers[cellPopulationMarkersRow, "CD8_positive"]
        )[, metaClustername], markerName] <-
          cellPopulationMarkers[cellPopulationMarkersRow, "name"]
      }
    }

    results <-
      merge(results, markerResults[, c(metaClustername, markerName)], by = metaClustername)
  }
  fwrite(
    results,
    paste0(
      "./data/",
      directoryName,
      '/clusteringOutput/',
      'meta_clusters_flowsom_Stability.csv'
    ),
    row.names =  FALSE
  )
}

stringBreak <- function(string,
                        sep,
                        buffer=80,
                        add=FALSE,
                        accomHash = FALSE,
                        breakchr="\n"
){

  if(buffer>90){ #Warn when output strings will certainly be too long for Mplus
    warning(paste0("All lines in the Mplus .inp file must contain less than 90 characters. ",
                   "Decrease the value of the buffer argument to avoid issues with Mplus parsing strings."))
  }

  if(buffer>89 && add==TRUE){ #Warn when output strings will certainly be too long for Mplus
    warning(paste0("All lines in the Mplus .inp file must contain less than 90 characters.",
                   "Appending the semicolon at the end of this strings may exceed this limit when buffer=90",
                   "Decrease the value of the buffer argument to avoid issues with Mplus parsing strings."))
  }

  if(!is(string,"character")){ #Return error if supplied an object not of character class
    stop("The object supplied to the string argument does not have character class, please convert to character and try again.")
  }

  #Create pattern for string subsetting, broadly: Match to last instance of `sep` within string
  pattern <- paste0(sep,"[^",sep,"]+$")


  #Check to see if the string exceeds the allowed buffer, and adjust if so
  if(nchar(string)>=buffer){

    #if the string is shorter than 2* buffer, subset only once
    if(nchar(string)<buffer*2){

      stringstart<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(stringstart,"", string)                   #Then, drop the stringstart portion from string, which may be the end of the string

    } else if(nchar(string)>=buffer*2){     #Subset using a while loop if the string >= buffer*2

      stringlist<- vector(mode="list")      #create empty list to store string subsets
      while(nchar(string)>=buffer){         #loop until stringend no longer exceeds buffer length
        subs<- sub(pattern,"", substr(string,1,buffer))              #Subset string to the first 1:buffer characters and break as specified in pattern
        stringlist[[length(stringlist)+1]] <- subs                   #Append subs object as next element of stringlist
        string <- sub(subs,"", string)                               #Drop the current subset portion from string
      } # end loop

      stringstart <- paste0(unlist(stringlist),collapse=breakchr) #paste unlisted stringlist and final string portion, collapsing with line breaks

    }

    #If there is not enough room for addition of a 32-character hash to the string, subset string once more
    if(accomHash==TRUE && buffer-nchar(string)<=33){
      pattern <- paste0(sep,"[^",sep,"]+$")

      fitHash<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(fitHash,"", string)                   #Then, drop the fitHash portion from string

      string <- paste0(fitHash,breakchr,string)
    }

    stringout<- paste0(stringstart,breakchr,string) #Paste string with linebreak between start and end


  } else {
    # #If there is not enough room for addition of a 32-character hash to the string, subset string
    if(accomHash==TRUE && buffer-nchar(string)<=33){
      fitHash<- sub(pattern,"", substr(string,1,buffer))  #Subset string to the first 1:buffer characters and break this substring as specified in pattern
      string <- sub(fitHash,"", string)                   #Then, drop the fitHash portion from string

      stringout <- paste0(fitHash,breakchr,string)
      message("The string was adjusted to accomodate the 32-character MD5 hash")
    } else {

      #If the string characters are already < buffer (taking into account accomHash), do nothing and return message
      message("The string already contains fewer characters than the buffer value defined, no changes have been made")
      stringout <- string

    }
  }

  #If add=TRUE, append a semicolon at the end of stringout, indicating the end of Mplus line
  if(add==TRUE){
    stringout <- paste0(stringout,";")
  }

  return(stringout)

}

reduceBoostrappedDataSize <- function(directoryName, iterations) {
  df <- fread(file=paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'))
  df <- as.data.frame(df)

  results <- data.frame(clusters_flowsom = unique(df[, "clusters_flowsom"]))

  for (number in seq(iterations)) {
    message(paste0("\nIteration: ", number))

    metaClustername <- paste0("meta_clusters_flowsom", number)
    clusterName <- paste0("clusters_flowsom", number)
    results[, metaClustername] <- NA

    subselectedDf <- na.omit(df[, c(clusterName, metaClustername)])

    for(cluster in unique(subselectedDf[, clusterName])){
      message(paste0("Cluster: ", cluster))
      metaCluster <- unique(subselectedDf[subselectedDf[, clusterName] == cluster, metaClustername])
      results[results[, "clusters_flowsom"] == cluster, metaClustername] <- metaCluster
    }

    df <- df[, colnames(df) != metaClustername]
  }

  fwrite(results, paste0("./data/", directoryName, '/clusteringOutput/', 'meta_clusters_flowsom_Stability.csv'), row.names =  FALSE)
  fwrite(df, paste0("./data/", directoryName, '/clusteringOutput/flowsomClusterStability.csv'), row.names =  FALSE)
}
