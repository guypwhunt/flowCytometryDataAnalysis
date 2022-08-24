###############################################################################################
##' Plots of CyTOF Analysis Ozlem
##' 30-01-2021 Johannes Schroth
###############################################################################################

## Setup ----
#rm(list = ls())

setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')


library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(cowplot)

## Heatmap ----
meta <- read.csv('Datatables/metaclusters.csv')
heat <- df

heat[heat$metaclusters == 'Excluded','clusters'] <- 'NA'

heat_mat <- heat %>%
  dplyr::select(panel$antigen[panel$marker_class == 'type'], clusters) %>%
  group_by(clusters) %>%
  summarise_all(median) %>%
  remove_rownames() %>%
  column_to_rownames('clusters') %>%
  as.matrix() %>%
  pheatmap:::scale_rows()

metacluster_colours <- as.vector(colortools::wheel('#ce302e', num = length(unique(heat$metaclusters))))
names(metacluster_colours) <- unique(heat$metaclusters)
metacluster_colours['Excluded'] <- 'grey50'

row_order <- meta$celltype[match(rownames(heat_mat), meta$clusters)]
row_order[is.na(row_order)] <- 'Excluded'
left_anno <- HeatmapAnnotation(Celltype = row_order, which = 'row', col=list('Celltype' = metacluster_colours), show_legend = F, show_annotation_name = F)
right_anno <- rowAnnotation(text1 = anno_text(paste0(round(matrix(table(heat$clusters)/nrow(heat)*100), 1),'%')))

p1 <- heat_mat %>%
  Heatmap(rect_gp = gpar(col = "white", lwd = 2),column_title_rot = 90,
          column_names_rot = 45, show_heatmap_legend = T,
          row_names_side = 'left', heatmap_legend_param = list(title = 'Z-score'),
          left_annotation = left_anno,
          row_split = row_order, row_gap = unit(2, "mm"), border = T, 
          row_title_rot = 0,
          col = circlize::colorRamp2(seq(-4, 4, length = 3), c('#4575B4', 'white', '#D73027')),
          right_annotation = right_anno)


png('Figures/Raw Plots/Heatmap.png', width = 4000, height = 2000, res = 320)
p1
dev.off()

## tSNE plots ----

clus_ds <- ds_group %>%
  dplyr::group_by(metaclusters) %>%
  dplyr::summarise(tSNE1 = median(tSNE1), 
                   tSNE2 = median(tSNE2))

p2 <- ds_group %>%
ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = as.factor(metaclusters)), alpha = 0.5, colour="black", pch = 21) +
  geom_label_repel(data = clus_ds, aes(label = metaclusters), alpha=0.9, colour = 'black', box.padding = 0.5, fontface = 'bold',segment.curvature = -0.1, segment.ncp = 3, segment.angle = 20) +
  theme_classic() +
  theme(legend.position = 'none', axis.text=element_blank(), axis.ticks=element_blank(), axis.title = element_text(face="bold", size = 14),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1)) +
  scale_fill_manual(values = metacluster_colours)


png('Figures/Raw Plots/tSNE.png', width = 2500, height = 2500, res = 320)
p2
dev.off()

### T cell clustering ----

p5 <- ggplot() +
  geom_point(data = ds_group[!ds_group$metaclusters %in% grep('T', ds_group$metaclusters, value = T),], aes(tSNE1, tSNE2), alpha = 0.5, colour="black", pch = 21, fill = 'grey50') +
  geom_point(data = ds_group[ds_group$metaclusters %in% grep('T', ds_group$metaclusters, value = T),], aes(tSNE1, tSNE2), alpha = 0.5, colour="black", pch = 21, fill = '#FF7878') +
  theme_classic() +
  theme(legend.position = 'none', axis.text=element_blank(), axis.ticks=element_blank(), axis.title = element_text(face="bold", size = 14),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1))


ggsave('Figures/Raw Plots/Tcells_highlighted_tSNE.png', p5, width = 10, height = 10, dpi = 500)

clus_ds <- tcell_ds_group %>%
  dplyr::group_by(metaclusters) %>%
  dplyr::summarise(tSNE1 = median(tSNE1), 
                   tSNE2 = median(tSNE2))

unique(clus_ds$metaclusters)
tcell_colours <- c("#C4D5F7","#9AABD5","#7081B3","#1B2C6F", "#F9D499","#EFAE6B","#E4873E","#DA6110")

names(tcell_colours) <- c('CD8 T N', 'CD8 T CM', 'CD8 T EFF', 'CD8 T SEN',
                          'CD4 T N', 'CD4 T CM', 'CD4 T EFF', 'CD4 T SEN')

p6 <- ggplot(tcell_ds_group, aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = as.factor(metaclusters)), alpha = 0.4, colour="black", pch = 21) +
  geom_density_2d(data = tcell_ds_group, aes(tSNE1, tSNE2), colour = 'black', alpha = 0.6) +
  scale_fill_manual(values = metacluster_colours) +
  facet_grid(~ factor(Group, c('CTRL', 'ALS')), labeller = labeller(Group = c('CTRL' = 'Control', 'ALS' = 'ALS'))) +
  ggpubr::theme_pubclean() +
  geom_label_repel(data = clus_ds, aes(label = metaclusters), alpha=0.9, colour = 'black', box.padding = 0.5, fontface = 'bold',segment.curvature = 0.1, segment.ncp = 10, segment.angle = 20) +
  theme(legend.position = 'none', axis.text=element_blank(), axis.ticks=element_blank(), axis.title = element_text(face="bold", size = 14),
        axis.line.y = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.4, "cm"), ends = "last"), size = 1),
        strip.text.x = element_text(size = 12, face = "bold"), strip.text.y = element_text(size = 12, face = "bold")) +
  expand_limits(x = 60, y = 20)
p6

ggsave('Figures/Raw Plots/Tcell_tSNE.png', p6, width = 16, height = 8, dpi = 500)

counts <- counts %>%
  mutate_all(function(x)x/sum(x)*100)

p4 <- counts %>%
  as.data.frame() %>%
  rownames_to_column('Clusters') %>%
  gather(key = 'sample_id', value = 'freq', -Clusters) %>%
  mutate(Group = md[match(.$sample_id, md$sample_id), 'Group']) %>%
  mutate(Group, Group = plyr::revalue(Group, c('CTRL' = 'HC'))) %>%
  mutate(Group = factor(Group, c('HC', 'ALS'))) %>%
  filter(Clusters == 'CD4 T SEN') %>%
  ggplot(aes(Group, freq)) +
  geom_boxplot(fill = metacluster_colours[names(metacluster_colours) == 'CD4 T SEN'], outlier.alpha = 0, alpha = 0.5) +
  geom_jitter(position = position_jitter(0.1)) +
  theme_pubclean() +
  facet_grid(~Clusters) +
  ylab('% live cells') +
  geom_bracket(xmin = "HC", xmax = "ALS", y.position = 10.5, inherit.aes = F,
               label = "p = 0.0017 (FDR-adjusted p < 0.05)", size = 0.5, label.size = 4, 
               tip.length = c(0.02, 0.02)) +
  theme(legend.position = 'none', 
        text = element_text(size=20),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), 
        strip.background = element_rect(fill = "grey90", colour = "black", size = 0.5)) +
  scale_y_continuous(breaks = seq(0,20,2.5), expand = c(0.01,1))

p4
ggsave('Figures/Raw Plots/Percent CD4 T senescent v2.png', p4, width = 6, height = 6, dpi = 500)




























