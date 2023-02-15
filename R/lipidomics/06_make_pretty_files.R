library(limma)
library(umap)
library(dplyr)
library(ggplot2)

figureDirectory <- "data/lipidomics/Results/"

fileNames <- list.files(figureDirectory, pattern =".csv")

fileNamesPath <- paste0(figureDirectory, fileNames)

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

unique(combinedDf$experiment)

combinedDf <- filterLipidomicsComparisons(combinedDf)

colnames(combinedDf)

combinedDf <- combinedDf[, c(
  "X",
  "experiment",
  "logFC",
  "AveExpr",
  "t",
  "B",
  "P.Value",
  "adj.P.Val"
)]

combinedDf$logPValue = 0-log10(combinedDf$P.Value)
combinedDf$logFDRPValue = 0-log10(combinedDf$adj.P.Val)


colnames(combinedDf) <- c(
  "Lipid",
  "Comparison",
  "Log Fold Change",
  "Average Expression",
  "T Statistic",
  "B Statistic (Log Odds)",
  "Raw P-Value",
  "FDR Adjusted P-Value",
  "Minus Log (Raw P-Value)",
  "Minus Log (FDR Adjusted P-Value)"
)

dir.create("data/pValueAdjustmentsResults", showWarnings = FALSE)

fwrite(combinedDf,
       paste0("data/pValueAdjustmentsResults/lipidomics.csv"))
