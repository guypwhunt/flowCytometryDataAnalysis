library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

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
  df$fdr_label <- gsub(' Negative/Low', '-', df$fdr_label)

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
    xlim(-0.5, 0.5) +
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
  df$fdr_label <- gsub(' Negative/Low', '-', df$fdr_label)

  #df[df$fdr_adjusted_p_val > 0.05, "fdr_label"] <- "Not Significant"

  df$minus_log_p_val <- 0 - log(df$p_val)
  if (exists("combinedDf")) {
    combinedDf <- rbind(combinedDf, df)
  } else {
    combinedDf <- df
  }
}

unique(combinedDf$fdr_label)

mycolors <- list(
  "HLA-DR- Activated CD11b+ Classical Monocytes" = "#f94144",
  "HLA-DR- Classical Monocytes" = "#f3722c",
  "HLA-DR- Intermediate Monocytes" = "#90be6d",
  "HLA-DR- Activated CD11b+ Intermediate Monocytes" = "#43aa8b",
  "HLA-DR- Non-Classical Monocytes" = "#577590",
  "HLA-DR- Activated CD11b+ Non-Classical Monocytes" = "#277da1",
  "Not Significant" = "black"
  )

labelValues <- c("ALS vs Healthy Controls",
                 "Fast Progressors vs Healthy Controls",
                 "Fast Progressors vs Slow Progressors",
                 "Slow Progressors vs Healthy Controls")

combinedDf$fdr_label <- as.factor(combinedDf$fdr_label)
combinedDf$comparison <- as.factor(combinedDf$comparison)

p <- ggplot(data = combinedDf,
            aes(
              x = logFC,
              y = minus_log_p_val,
              col = fdr_label,
              label = fdr_label
            )) +
  geom_point(
    aes(shape=comparison)
    ) +
  theme_minimal() +
  #scale_colour_viridis_d() +
  scale_colour_manual(values = mycolors) +
  xlab("log2(Fold Change)") + ylab("-log10(P-value)") +
  xlim(-0.5, 0.5) +
  ylim(0, 10) +
  scale_shape_discrete(labels = labelValues) +
  guides(colour = guide_legend("Cell Types"),
         shape = guide_legend("Comparison")) +
  theme(axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
  #geom_text_repel()
print(p)

filteredcombinedDf <- combinedDf

#filteredcombinedDf <- filteredcombinedDf[filteredcombinedDf$fdr_adjusted_p_val <= 0.05,]
filteredcombinedDf <- filteredcombinedDf[filteredcombinedDf$fdr_label != "Monocytes",]


filteredcombinedDf$fdr_label <- gsub(" Classical ", " \n Clasical ", filteredcombinedDf$fdr_label)
filteredcombinedDf$fdr_label <- gsub(" Non-Classical ", " \n Non-Classical ", filteredcombinedDf$fdr_label)
filteredcombinedDf$fdr_label <- gsub(" Intermediate ", " \n Intermediate ", filteredcombinedDf$fdr_label)
filteredcombinedDf$fdr_label <- paste0(filteredcombinedDf$fdr_label, " (", filteredcombinedDf[, "cluster_id"], ")")

filteredcombinedDf[, c("cluster_id", "fdr_label")]

filteredcombinedDf$comparison <- str_replace_all(filteredcombinedDf$comparison, "_", " ")


head(filteredcombinedDf)

ggplot(filteredcombinedDf, aes(x = fdr_label, y = minus_log_p_val,
                      color = as.factor(comparison), size = logFC)) +
  geom_point() +
  xlab("Monocyte Populations") + ylab("-log10(P-value)") +
  guides(colour = guide_legend("Comparisons",
                               override.aes = list(size=5)),
         size = guide_legend("Log2(Fold Change)")) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  scale_colour_viridis_d(alpha = 0.75) +
  ylim(0, 10) +
  geom_hline(yintercept = 0-log10(0.00001))



setwd(workingDirectory)
