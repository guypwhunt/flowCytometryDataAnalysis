try(source("R/01_functions.R"))

loadlibraries()

workingDirectory <- getwd()

directoryNames <- c("tCells_ISO", "tCells_FMO")

columnNames <-
  c(
    "CD127.BV510.A",
    "CD3.BV605.A",
    "CD8.BV650.A",
    "CD25.BV786.A",
    "FoxP3.PE.A",
    "CD45RO.PE.CF595.A",
    "CD4.PerCP.Cy5.5.A",
    "GPR32.AF488.A",
    "Zombie.NIR.A"
  )

automatedcofactors <- c("CD127", "CD3","CD8",
                        "CD25", "FoxP3","CD45RO",
                        "CD4", "GPR32", "Zombie")

for (directoryName in directoryNames) {
  setwd(paste0("./data/isotypes/", directoryName))

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
      print(head(raw[IND,]))
      cat("----\n")
    }
    return(unique(raw))
  }

  # Read all:
  dfs <- sapply(filenames, read.flow_csv, simplify = FALSE)

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

  # rewrite to flowset
  dfs_fs <- as(dfs_ff, "flowSet")

  #auto
  dfs_fs_t_auto <- transFlowVS(dfs_fs, channels = colnames(dfs_fs),
                               cofactor = rep(150, length(colnames(dfs_fs))))

  n <- 1
  saveFlowSetAsCsv <- function(df, columnNames, filenames) {
    df <- as.data.frame(exprs(df))
    write.csv(df[,columnNames], paste0("./transformedData/", filenames[n]), row.names = FALSE)
    n <<- n +1
  }

  fsApply(dfs_fs_t_auto, saveFlowSetAsCsv, columnNames, filenames)

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
