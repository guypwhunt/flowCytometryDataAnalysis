try(source("R/01_functions.R"))

loadlibraries()
library(gtools)


workingDirectory <- getwd()

results <- data.frame(
  Type = character(),
  FirstFile = character(),
  SecondFile = character(),
  Marker = character(),
  P_value = double(),
  Fold_Change = double(),
  Minus_P_value = double(),
  Adjusted_P_value = double(),
  Minus_Adjusted_P_value = double()
)

directoryNames <-
  c("bCells_ISO","monocytes_ISO", "tCells_ISO", "senescence_ISO")

columnNames <-
  c("GPR18", "Chem23")

for (directoryName in directoryNames) {
  setwd(paste0("./data/isotypesGPR18/", directoryName, "/transformedData"))

  filenames <- list.files(pattern = ".csv")

  lengthfilenames <- length(filenames)

  columnNamesDf <- read.csv(filenames[1])

  for (columnName in colnames(columnNamesDf)) {
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

      xdf <- data.frame(group,marker)
      firstPair <- xdf[xdf$group == filenames[i-2], "marker"]
      secondPair <- xdf[xdf$group == filenames[i-1], "marker"]
      foldChange <- foldchange(median(firstPair), median(secondPair))

      results[nrow(results) + 1,] <-
        c(directoryName, filenames[lownumber], filenames[highNumber], columnName, test[1], foldChange, 0,0,0)
    }

  }
  setwd(workingDirectory)


  results[,"P_value"] <- as.double(results[,"P_value"])
  results[results[,"P_value"] == 0, "P_value"] <- 5e-324
  results[, "Minus_P_value"] <- 0-log10(results[,"P_value"])
  results["Adjusted_P_value"] <- p.adjust(results[,"P_value"], method = "fdr", n = nrow(results))
  results["Minus_Adjusted_P_value"] <- 0-log10(results["Adjusted_P_value"])
  dir.create("./data/isotypesGPR18/pairwiseTTestResults/", showWarnings = FALSE)
  write.csv(results, paste0("./data/isotypesGPR18/pairwiseTTestResults/", directoryName, ".csv"))
}
