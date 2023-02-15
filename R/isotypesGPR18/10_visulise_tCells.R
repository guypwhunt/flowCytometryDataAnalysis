try(source("R/01_functions.R"))

loadlibraries()
library(gtools)

workingDirectory <- getwd()

directoryName <- "tCells_ISO"

setwd(paste0("./data/isotypes/", directoryName))

# Create an 'output' folder
gc()
dir.create("results", showWarnings = FALSE)
dir.create("transformedData", showWarnings = FALSE)

gc()

# Find file names of .csv files in the current working directory:
filenames <- list.files(pattern = ".csv")
isoFilenames <- list.files(pattern = "ISO.csv")
filenames <- filenames[!filenames %in% isoFilenames]

filenames <- filenames[1:4]
isoFilenames <- isoFilenames[1:4]

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

outliers_removed_dfs <- list()

i <- 1
for (df in dfs) {
  xx <- df[,c("CD127.BV510.A", "CD25.BV786.A")]

  cd11bQuantiles <- quantile(xx$CD127.BV510.A, c(0.001, 0.999))

  xx <- xx[xx$CD127.BV510.A > min(cd11bQuantiles) & xx$CD127.BV510.A < max(cd11bQuantiles), ]

  cd11bActivatedQuantiles <- quantile(xx$CD25.BV786.A , c(0.01, 0.99))

  xx <- xx[xx$CD25.BV786.A > min(cd11bActivatedQuantiles) & xx$CD25.BV786.A < max(cd11bActivatedQuantiles), ]

  xx <- list(xx)
  names(xx) <- names(dfs[i])

  i <- i+1
  outliers_removed_dfs <- append(outliers_removed_dfs, xx)
}


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

dfs_ff_outliers <- sapply(outliers_removed_dfs, function(x)
  csv_2_ff(x), simplify = FALSE)

dfs_outliers_fs <- as(dfs_ff_outliers, "flowSet")

# rewrite to flowset
dfs_fs <- as(dfs_ff, "flowSet")

colnames(dfs_fs)

dfs_fs_t_auto <- transFlowVS(dfs_fs, channels = c("CD127.BV510.A", "CD25.BV786.A"),
                             cofactor = c(500, 150))

figureDirectory <- paste0(getwd(), "/results/")

library(ggcyto)
ggcyto(dfs_fs, aes_string("CD127.BV510.A", "CD25.BV786.A")) +
  geom_hex(bins = 128)
ggcyto(dfs_outliers_fs, aes_string("CD127.BV510.A", "CD25.BV786.A")) +
  geom_hex(bins = 128)
ggcyto(dfs_fs_t_auto, aes_string("CD127.BV510.A", "CD25.BV786.A")) +
  geom_hex(bins = 128)

columnNameFormula <- as.formula(paste(" ~ ", "CD127.BV510.A"))
densityplot(columnNameFormula, dfs_fs, main = "auto")
densityplot(columnNameFormula, dfs_outliers_fs, main = "auto")
densityplot(columnNameFormula, dfs_fs_t_auto, main = "auto")

columnNameFormula <- as.formula(paste(" ~ ", "CD25.BV786.A"))
densityplot(columnNameFormula, dfs_fs, main = "auto")
densityplot(columnNameFormula, dfs_outliers_fs, main = "auto")
densityplot(columnNameFormula, dfs_fs_t_auto, main = "auto")

