try(source("R/01_functions.R"))

loadlibraries()

workingDirectory <- getwd()

directoryNames <- c("bCells_ISO", "monocytes_ISO", "senescence_ISO", "tCells_ISO")
# directoryName <- "tCells_ISO"

columnNames <-
  c("GPR32-AF488-A", "GPR32 - AF488-A", "GPR18 - AF488-A", "receptor - AF488-A", "Chem23 - AF647-A", "FPRL1 - AF647-A", "Chem23 - APC-A", "receptor - AF647-A", "FPRL1-AF647-A")

prettyColumnNames <- c("GPR18", "GPR18", "GPR18", "GPR18", "Chem23", "Chem23", "Chem23", "Chem23", "Chem23")

print(as.data.frame(columnNames, prettyColumnNames))

automatedcofactors <- c(2.166595e+01, 10.386289)

for (directoryName in directoryNames) {
  setwd(paste0("./data/isotypesGPR18/", directoryName))

  # Create an 'output' folder
  gc()
  dir.create("results", showWarnings = FALSE)
  dir.create("transformedData", showWarnings = FALSE)

  gc()

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  ## Defining a function to read a flow cytrometry file in csv format:
  # Each row is a cell, each column is a parameter. In our experience, the
  # flow cytometers sometimes output duplicate entries (listing the same cell
  # twice), we remove these and report.
  # Please check how your csv file is separated and adjust the sep argument
  # in the function if necessary. In this example we import a semicolon
  # separated file.
  read.flow_csv <- function(pathIN) {
    message(pathIN)
    raw <-
      fread(
        pathIN#,
        # ,
        # sep = ",",
        # header = TRUE,
        # stringsAsFactors = FALSE,
        # nrows = 100
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
    # if ("x-A" %in% colnames(df)) {
    #   print(colnames(df))
    # }
    # print("")

    xlist <- append(xlist, colnames(df))
    xcolumnNames <- columnNames[columnNames %in% colnames(df)]
    xprettyColumnNames <- prettyColumnNames[columnNames %in% colnames(df)]

    # print(columnNames[!columnNames %in% colnames(df)])

    df <- df[, xcolumnNames]
    colnames(df) <- xprettyColumnNames
    dfs[[x]]  <- df

    # if(length(xcolumnNames) < 9){
    #   print(names(dfs[x]))
    #   print(xcolumnNames)
    # }
  }

  unique(xlist)

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
  dfs_ff = sapply(dfs, function(x)
    csv_2_ff(x), simplify = FALSE)

  # rewrite to flowset
  dfs_fs <- as(dfs_ff, "flowSet")

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels = colnames(dfs_fs),
                               cofactor = automatedcofactors)

  n <- 1
  saveFlowSetAsCsv <- function(df, columnNames, filenames) {
    df <- as.data.frame(exprs(df))
    write.csv(df[,columnNames], paste0("./transformedData/", filenames[n]), row.names = FALSE)
    n <<- n +1
  }

  fsApply(dfs_fs_t_auto, saveFlowSetAsCsv, unique(prettyColumnNames), filenames)

  # Pre-Normalized Plots
  flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

  figureDirectory <- paste0(getwd(), "/results/")

  for (columnName in colnames(dfs_fs)) {
    columnNameFormula <- as.formula(paste(" ~ ", columnName))

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
  }


  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}
