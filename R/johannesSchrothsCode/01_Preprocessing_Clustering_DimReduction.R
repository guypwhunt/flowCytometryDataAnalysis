###############################################################################################
##' CyTOF data analysis of ALS patient PBMCs
##' 29-03-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------------------------------------------------

#rm(list = ls())

set.seed(1234)

setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

library(flowCore)
library(tidyverse)
library(cytofkit2)
library(rsvd)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(survival)
library(survminer)
library(Hmisc)
library(ggcorrplot)
library(emmeans)
library(lme4)
library(multcomp)
library(FastPG)

source('~/development/FIt-SNE-master/fast_tsne.R', chdir=T)
reticulate::use_python('/usr/local/bin/python3.9', required = T)

## Data Import and Preprocessing --------------------------------------------------------------------------------------------------------

md <- read.csv('./Datatables/metadata.csv')

panel <- read.csv('Datatables/panel.csv')

fs <- read.flowSet(files = as.vector(md$file_name),
                   transformation = F,
                   truncate_max_range = F)

sampleNames(fs) <- sub('.*(BLT).*(\\d{5}).*', '\\1\\2' , sampleNames(fs))
colnames(fs) <- sub('(.*)..v.', '\\1', parameters(fs[[1]])$desc)
fs <- fs[, match(panel$fcs_colname, colnames(fs))]
colnames(fs) <- panel$antigen[match(panel$fcs_colname, colnames(fs))]

fs_transformed <- fsApply(fs, function(x, cofactor = 5) {
  expr <- exprs(x)
  expr <- asinh(expr / cofactor)
  exprs(x) <- expr
  x
})

sample_length <- fsApply(fs_transformed, nrow)

df <- data.frame(fsApply(fs_transformed, flowCore::exprs),
                 'sample_id' = rep(rownames(sample_length), sample_length))

df <- df %>%
  mutate(Group = md[match(df$sample_id, md$sample_id), 'Group'])

#rm(fs, fs_transformed, fs_corrected)

## Clustering ---------------------------------------------------------------------------------------------------------------------------

df['clusters'] <- df %>%
  select(panel$antigen[panel$marker_class == 'type']) %>%
  as.matrix() %>%
  fastCluster(k = 30, num_threads = 4) %>%
  .$communities

# Exclude clusters of <1%
table(df$clusters)/nrow(df)*100
sum(table(df$clusters)[table(df$clusters)/nrow(df) < 0.01]/nrow(df))

df %>%
  dplyr::select(panel$antigen[panel$marker_class == 'type'], clusters) %>%
  group_by(clusters) %>%
  summarise_all(median) %>%
  dplyr::filter(clusters %in% as.numeric(names(which(table(df$clusters)/nrow(df) > 0.01)))) %>%
  remove_rownames() %>%
  column_to_rownames('clusters') %>%
  as.matrix() %>%
  pheatmap:::scale_rows() %>%
  Heatmap(rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 45, show_heatmap_legend = T,
          row_names_side = 'left', heatmap_legend_param = list(title = 'Z-score'),
          col = circlize::colorRamp2(seq(-4, 4, length = 3), c('#4575B4', 'white', '#D73027')))

meta <- read.csv('./Datatables/metaclusters.csv')

df <- df %>%
  mutate(metaclusters = meta[match(clusters, meta$clusters), 'celltype'])

## Dimensionality Reduction -------------------------------------------------------------------------------------------------------------

ds_group <- df %>%
  filter(metaclusters != 'Excluded') %>%
  groupdata2::balance(., 50000, 'Group')

ds_group[c('tSNE1', 'tSNE2')] <- ds_group %>%
  select(panel$antigen[panel$marker_class == 'type']) %>%
  as.matrix() %>%
  fftRtsne(perplexity = nrow(.)/100,
           learning_rate = nrow(.)/12,
           theta = 0.5, 
           max_iter = 1000,
           rand_seed = 1234)

## T cell extraction, reclustering and analysis -----------------------------------------------------------------------------------------

tcell_df <- df %>%
  filter(metaclusters %in% grep('CD4|CD8', unique(df$metaclusters), value = T))

tcell_ds_group <- tcell_df %>%
  filter(metaclusters != 'Excluded') %>%
  groupdata2::balance(., 50000, 'Group')

tcell_ds_group[c('tSNE1', 'tSNE2')] <- tcell_ds_group %>%
  dplyr::select(panel$antigen[panel$marker_class == 'type']) %>%
  as.matrix() %>%
  fftRtsne(perplexity = nrow(.)/100,
           learning_rate = nrow(.)/12,
           theta = 0.5, 
           max_iter = 1000,
           rand_seed = 1234)

tcell_ds_group %>%
  ggplot(aes(tSNE1, tSNE2, colour = metaclusters)) +
  geom_point()

save.image(file = './R/CyTOF_Data.RData')


# Get the median marker expression per sample
library(dplyr)


expr_median_sample_tbl <- df %>%
  filter(metaclusters == 'CD4 T SEN') %>%
  group_by(sample_id) %>%
  summarise_at(panel$antigen[panel$marker_class == 'type'], median)

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)
mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)

ggdf$condition <- md$Group[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, fill = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw()


