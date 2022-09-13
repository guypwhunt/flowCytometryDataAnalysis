library(dplyr)
library(ggplot2)
library(ggrepel)

workingDirectory <- getwd()

setwd(paste0("./data/pValueAdjustmentsResults/"))

files <- c(
  "meta_clusters_flowsomcaseControlVisits1AllCellsDifferentialStatesStatisticscsvClusters.csv",
  "meta_clusters_flowsomcaseControlVisits1FastAllCellsDifferentialStatesStatisticscsvClusters.csv",
  "meta_clusters_flowsomcaseControlVisits1SlowAllCellsDifferentialStatesStatisticscsvClusters.csv",
  "meta_clusters_flowsomfastSlowVisits1AllCellsDifferentialStatesStatisticscsvClusters.csv"
)

for (file in files) {
  df <-
    read.csv(file)

  df <- df[df$panel == "monocytes", ]

  df$fdr_label <- df$typeOfCells
  #df$fdr_label <- gsub('Monocytes', '', df$fdr_label)
  df$fdr_label <- gsub(' Negative/Low', '-', df$fdr_label)
  #paste0(df$typeOfCells, " (Cluster ", df$cluster_id, ")")
  #df$fdr_label[1] <- "Monocytes"

  df[df$fdr_adjusted_p_val > 0.05, "fdr_label"] <- NA

  df$minus_log_p_val <- 0 - log(df$p_val)

  mycolors <- data.frame(DOWN = "blue",
                         UP = "red",
                         NO = "black")

  p <- ggplot(data = df,
              aes(
                x = logFC,
                y = minus_log_p_val,
                col = fdr_diff_expressed,
                label = fdr_label
              )) +
    geom_point() +
    theme_minimal() +
    scale_colour_manual(values = mycolors) +
    xlab("log2(Fold Change)") + ylab("-log10(P-value)") +
    xlim(-0.05, 0.05) +
    ylim(0, 10) +
    geom_text_repel(size = 2, force_pull = 0, max.time = 2) +
    theme(legend.position = "none")
  print(p)

}

comparisons <- c(
  "ALS_vs_Healthy_Controls",
  "Fast_Progressors_vs_Healthy_Controls",
  "Slow_Progressors_vs_Healthy_Controls",
  "Fast_Progressors_vs_Slow_Progressors"
)

filesDf <- data.frame(files, comparisons)

for (file in seq(nrow(filesDf))) {
  df <-
    read.csv(filesDf[file, "files"])

  df <- df[df$panel == "monocytes", ]
  df$comparison <- filesDf[file, "comparisons"]


  df$fdr_label <- df$typeOfCells
  #df$fdr_label <- gsub('Monocytes', '', df$fdr_label)
  df$fdr_label <- gsub(' Negative/Low', '-', df$fdr_label)
  #paste0(df$typeOfCells, " (Cluster ", df$cluster_id, ")")
  #df$fdr_label[1] <- "Monocytes"

  df[df$fdr_adjusted_p_val > 0.05, "fdr_label"] <- NA
  df[df$fdr_adjusted_p_val > 0.05, "comparison"] <- "Not_Significant"

  df$minus_log_p_val <- 0 - log(df$p_val)
  if (exists("combinedDf")) {
    combinedDf <- rbind(combinedDf, df)
  } else {
    combinedDf <- df
  }
}

mycolors <- data.frame(
  "ALS_vs_Healthy_Controls" = "red",
  "Slow_Progressors_vs_Healthy_Controls" = "yellow",
  "Fast_Progressors_vs_Healthy_Controls" = "green",
  "Fast_Progressors_vs_Slow_Progressors" = "blue",
  "Not_Significant" = "black"
  )

labelValues <- c("ALS vs Healthy Controls",
                 "Fast Progressors vs Healthy Controls",
                 #"Slow Progressors vs Healthy Controls" = "8ac926",
                 "Fast Progressors vs Slow Progressors" ,
                 "Not Significant")

combinedDf$fdr_label <- as.factor(combinedDf$fdr_label)
combinedDf$comparison <- as.factor(combinedDf$comparison)

p <- ggplot(data = combinedDf,
            aes(
              x = logFC,
              y = minus_log_p_val,
              col = comparison,
              label = comparison
            )) +
  geom_point() +
  theme_minimal() +
  scale_colour_manual(values = mycolors) +
  xlab("log2(Fold Change)") + ylab("-log10(P-value)") +
  geom_label_repel(size = 2) +
  theme(legend.position = "none") +
  xlim(-0.05, 0.05) +
  ylim(0, 10)
print(p)

setwd(workingDirectory)
