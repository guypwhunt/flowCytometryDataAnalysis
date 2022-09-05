try(source("R/01_functions.R"))

loadlibraries()

workingDirectory <- getwd()

directoryNames <- c("bCells_ISO", "bCells_FMO", "monocytes_ISO", "monocytes_FMO")

columnNames <-
  c("GPR32...AF488.A", "FPRL1...AF647.A")

for (directoryName in directoryNames) {
  i <- 1
  fileI <- 1

  setwd(paste0("./data/isotypes/", directoryName))

  # Create an 'output' folder
  gc()
  dir.create("results", showWarnings = FALSE)
  gc()

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  dfs = lapply(filenames, read.csv)

  results <- data.frame(
    fileName = character(),
    marker = character(),
    min = double(),
    median = double(),
    max = double()
  )

  for (df in dfs) {
    for (columnName in columnNames) {
      results[i, ] <- c(filenames[fileI], columnName,
                        min(df[, columnName]),
                        median(df[, columnName]),
                        max(df[, columnName]))
      i <- i + 1
    }
    fileI <- fileI + 1
  }

  write.csv(results, paste0("results/", directoryName, "summaryStats.csv"))

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}


####

directoryNames <- c("tCells_ISO", "tCells_FMO")

columnNames <-
  c("GPR32.AF488.A",
    "FPRL1.AF647.A")

for (directoryName in directoryNames) {
  i <- 1
  fileI <- 1

  setwd(paste0("./data/isotypes/", directoryName))

  # Create an 'output' folder
  gc()
  dir.create("results", showWarnings = FALSE)
  gc()

  # Find file names of .csv files in the current working directory:
  filenames <- list.files(pattern = ".csv")

  dfs = lapply(filenames, read.csv)

  results <- data.frame(
    fileName = character(),
    marker = character(),
    min = double(),
    median = double(),
    max = double()
  )

  for (df in dfs) {
    for (columnName in columnNames) {
      results[i, ] <- c(filenames[fileI], columnName,
                        min(df[, columnName]),
                        median(df[, columnName]),
                        max(df[, columnName]))
      i <- i + 1
    }
    fileI <- fileI + 1
  }

  write.csv(results, paste0("results/", directoryName, "summaryStats.csv"))

  tryCatch({
    setwd(workingDirectory)
  },
  error = function(cond) {
    setwd("..")
    setwd("..")
  })
}
