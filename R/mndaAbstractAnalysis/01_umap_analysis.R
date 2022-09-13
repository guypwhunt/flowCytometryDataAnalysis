library(dplyr)
library(ggplot2)
library(ggrepel)

workingDirectory <- getwd()

setwd("./data/monocytes/clusteringOutput/")

df <- read.csv('umapDf.csv')

metaFlowSomcellPopulations <-
  read.csv('meta_clusters_flowsomCellPopulations.csv')
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

df$meta_flowsom_cell_population <- gsub(' Negative/Low', '-', df$meta_flowsom_cell_population)

df[, "cell_population"] <- paste0(df[, "meta_flowsom_cell_population"], " (Cluster ", df[, "meta_clusters_flowsom"], ")")

#visualize and label clusters on umap
gc()
label_meta_flowsom_umap <-
  df %>% group_by(meta_flowsom_cell_population) %>%
  select(umap_1, umap_2) %>% summarize_all(mean)

gc()
plot <-
  ggplot(df, aes(
    x = umap_1,
    y = umap_2,
    color = as.factor(meta_flowsom_cell_population)
  )) + geom_point(size = 0.1) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "None"
  ) +
  geom_label_repel(aes(label = meta_flowsom_cell_population),
                   data = label_meta_flowsom_umap,
                   size = 2, force_pull = 0, max.time = 2) +
  labs(y= "UMAP 2", x = "UMAP 1")

print(plot)

setwd(workingDirectory)
