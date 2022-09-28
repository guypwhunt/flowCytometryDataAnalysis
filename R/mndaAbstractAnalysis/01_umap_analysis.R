library(dplyr)
library(ggplot2)
library(ggrepel)

workingDirectory <- getwd()

setwd("./data/monocytes/automatedCofactors_clusteringOutput/")

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

df$meta_flowsom_cell_population <- paste0(df$meta_flowsom_cell_population, " (", df[, "meta_clusters_flowsom"], ")")

#df$meta_flowsom_cell_population <- gsub(" Classical ", " \n Clasical ", df$meta_flowsom_cell_population)
#df$meta_flowsom_cell_population <- gsub(" Non-Classical ", " \n Non-Classical ", df$meta_flowsom_cell_population)
#df$meta_flowsom_cell_population <- gsub(" Intermediate ", " \n Intermediate ", df$meta_flowsom_cell_population)

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
  labs(y= "UMAP 2", x = "UMAP 1") +
  guides(colour = guide_legend("Monocyte Populations",
                               override.aes = list(size=5)))

print(plot)

label_meta_flowsom_umap <-
  df %>% group_by(meta_clusters_flowsom) %>%
  select(umap_1, umap_2) %>% summarize_all(mean)


gc()
plot <-
  ggplot(df, aes(
    x = umap_1,
    y = umap_2,
    color = as.factor(meta_clusters_flowsom)
  )) + geom_point(size = 0.1) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "None"
  ) +
  geom_label_repel(aes(label = meta_clusters_flowsom),
                   data = label_meta_flowsom_umap,
                   size = 3, force_pull = 0, max.time = 2) +
  labs(y= "UMAP 2", x = "UMAP 1")

print(plot)

setwd(workingDirectory)

library(gtools)
x <- c(1,2,3,4,5)
y <- c(2,2,2,2,2)
x <- median(x)
y <- median(y)

foldchange(x, y)
