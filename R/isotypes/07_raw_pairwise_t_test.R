try(source("R/01_functions.R"))

loadlibraries()

workingDirectory <- getwd()

results <- data.frame(
  Type = character(),
  FirstFile = character(),
  SecondFile = character(),
  Marker = character(),
  P_value = double()
)

directoryNames <-
  c("bCells_ISO", "bCells_FMO", "monocytes_ISO", "monocytes_FMO")

columnNames <-
  c("GPR32...AF488.A", "FPRL1...AF647.A")

for (directoryName in directoryNames) {
  setwd(paste0("./data/isotypes/", directoryName))

  filenames <- list.files(pattern = ".csv")

  lengthfilenames <- length(filenames)

  for (columnName in columnNames) {
    for (x in (seq(lengthfilenames / 2))) {
      highNumber <- 2 * x
      lownumber <- 2 * x - 1
      dfs = lapply(filenames[lownumber:highNumber], read.csv)

      marker <- c()
      group <- c()

      i <- lownumber
      for (df in dfs) {
        marker <- append(marker, df[, columnName])
        group <- append(group, rep(filenames[i], nrow(df)))
        i <- i + 1
      }

      test <-
        pairwise.t.test(marker, as.factor(group), p.adjust.method = "none")$p.value

      results[nrow(results) + 1,] <-
        c(directoryName, filenames[lownumber], filenames[highNumber], columnName, test[1])
    }

  }
  setwd(workingDirectory)
}

directoryNames <- c("tCells_ISO", "tCells_FMO")

columnNames <-
  c("GPR32.AF488.A",
    "FPRL1.AF647.A")


for (directoryName in directoryNames) {
  setwd(paste0("./data/isotypes/", directoryName))

  filenames <- list.files(pattern = ".csv")

  lengthfilenames <- length(filenames)

  for (columnName in columnNames) {
    for (x in (seq(lengthfilenames / 2))) {
      highNumber <- 2 * x
      lownumber <- 2 * x - 1
      dfs = lapply(filenames[lownumber:highNumber], read.csv)

      marker <- c()
      group <- c()

      i <- lownumber
      for (df in dfs) {
        marker <- append(marker, df[, columnName])
        group <- append(group, rep(filenames[i], nrow(df)))
        i <- i + 1
      }

      test <-
        pairwise.t.test(marker, as.factor(group), p.adjust.method = "none")$p.value

      results[nrow(results) + 1,] <-
        c(directoryName, filenames[lownumber], filenames[highNumber], columnName, test[1])
    }

  }
  setwd(workingDirectory)
}

results["Adjusted_P_value"] <- p.adjust(results[,"P_value"], method = "fdr", n = nrow(results))

write.csv(results, "./data/isotypes/rawSummaryStats/rawSummaryStats.csv")
