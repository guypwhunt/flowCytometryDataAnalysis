###############################################################################################
##' CyTOF data analysis of ALS patient PBMCs
##' 29-03-2021 Johannes Schroth
##' Supplementary Figure comparing manual and phenograph cluster abundance
###############################################################################################

## Setup --------------------------------------------------------------------------------------

#rm(list = ls())

setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

library(tidyverse)
library(ggrepel)
library(RColorBrewer)


p1_tmp <- df %>%
  dplyr::select(panel$antigen[panel$marker_class == 'type'], sample_id) %>%
  gather(key, value, -sample_id) %>%
  ggplot(aes(x = value, colour = sample_id)) +
  geom_density() +
  facet_wrap(~key, scales = 'free', nrow = 3) +
  theme_pubclean() +
  theme(legend.position = 'none')

ggsave('Figures/Supplementary Plots/Marker expression by Sample.png', p1_tmp, width = 10, height = 6, dpi = 500)

p1_tmp <- counts %>%
  rownames_to_column('Metaclusters') %>%
  gather(key = 'sample_id', value = 'freq', -Metaclusters) %>%
  mutate(Group = factor(md[match(.$sample_id, md$sample_id), 'Group'], levels = c('CTRL', 'ALS'))) %>%
  mutate(Group, Group = plyr::revalue(Group, c('CTRL' = 'HC'))) %>%
  ggplot(aes(x = sample_id, y = freq * 100, fill = Metaclusters)) +
  geom_bar(stat = 'identity', colour = 'black') +
  theme_pubclean() +
  scale_fill_manual(values = metacluster_colours) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~ Group, scales = 'free', space = 'free') +
  theme(legend.position = 'right', 
        text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border = element_rect(colour = "black", fill = NA), 
        strip.background = element_rect(fill = "grey90", colour = "black", size = 0.5)) +
  ylab('% live cells')

p1_tmp

ggsave('Figures/Supplementary Plots/Metacluster Frequency by Sample v2.png', p1_tmp, width = 15, height = 10, dpi = 500)

norm_melt <- ds_group %>%
  dplyr::select(panel$antigen, clusters, tSNE1, tSNE2) %>%
  gather(variable, value, -clusters, -tSNE1, -tSNE2)

unique(norm_melt$variable)

p3 <- list()

for (i in unique(norm_melt$variable)) {
  p3[[i]] <- ggplot(norm_melt[norm_melt$variable == i,], aes(tSNE1, tSNE2, colour = value)) +
    geom_point(alpha = 0.5) +
    theme_void() +
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 8)) +
    ggtitle(paste(i)) +
    scale_colour_gradient2(low="#4575B4", mid= 'white' , high='#D73027', midpoint=mean(norm_melt[norm_melt$variable == i,'value']))
}
get_legend()

cowplot::save_plot('./Figures/Supplementary Plots/marker expression matrix.png', cowplot::plot_grid(plotlist = p3, nrow = 4), base_width = 35, base_height = 20, units = "cm")


## Supplementary Fig 1 ------------------------------------------------------------------------

manual_gating_df <- read.csv('Datatables/Manual Gating/manual_gating_frequency.csv', header = T)
pheno_gating_df <- read.csv('Datatables/Results/Cluster Frequency.csv')

head(manual_gating_df)
head(pheno_gating_df)

manual_gating_df <- manual_gating_df %>%
  column_to_rownames('X') %>%
  t() %>%
  as.data.frame()

rownames(manual_gating_df) <- sub('..', ' +', rownames(manual_gating_df), fixed = T)
rownames(manual_gating_df) <- sub('.', ' ', rownames(manual_gating_df), fixed = T)
rownames(manual_gating_df) <- sub('.', ' ', rownames(manual_gating_df), fixed = T)
rownames(manual_gating_df)[7] <- 'Non-classical Monocyte'
rownames(manual_gating_df)
head(manual_gating_df)


manual_gating_df <- manual_gating_df %>%
  rownames_to_column('celltype') %>%
  gather(key = 'id', value = 'percentage', -celltype) %>%
  mutate(Gating = 'manual')

pheno_gating_df <- pheno_gating_df %>%
  dplyr::rename(celltype = X) %>%
  gather(key = 'id', value = 'percentage', -celltype) %>%
  mutate(Gating = 'phenograph',
         percentage = percentage * 100)

gating_df <- rbind(manual_gating_df, pheno_gating_df)
gating_df$celltype <- reorder(gating_df$celltype, -gating_df$percentage, mean)

gating_df$cellgroup <- NA
gating_df$cellgroup[grep('Monocyte', gating_df$celltype)] <- 'Monocytes'
gating_df$cellgroup[grep('B cell', gating_df$celltype)] <- 'B cells'
gating_df$cellgroup[grep('NK', gating_df$celltype)] <- 'NK cells'
gating_df$cellgroup[grep('CD4|CD8|Treg|gd T', gating_df$celltype)] <- 'T cells'

p1_tmp <- gating_df %>%
  dplyr::filter(!celltype %in% c('CD4 T EM', 'CD8 T EM')) %>%
  ggplot(aes(celltype, percentage, fill = Gating)) +
  geom_boxplot(alpha = 0.8) +
  ggpubr::theme_pubclean() +
  scale_fill_discrete(breaks = c('manual', 'phenograph'), labels = c('Manual Gating', 'Phenograph Clustering')) +
  geom_hline(data=tibble(f=1, y=1*Inf), aes(yintercept=y), col="black", size = 1) +
  facet_grid(~cellgroup, scales = 'free', space = "free_x") +
  labs(x = NULL, y = '% Live Cells') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.background=element_rect(fill="white"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=20))

ggsave('Figures/Supplementary Plots/Manual and Phenograph Gating Comparison.png', p1_tmp, width = 15, height = 10, dpi = 500)




p1_tmp <- counts %>%
  as.data.frame() %>%
  rownames_to_column('Clusters') %>%
  gather(key = 'sample_id', value = 'freq', -Clusters) %>%
  mutate(Group = md[match(.$sample_id, md$sample_id), 'Onset']) %>%
  filter(Group %in% c('B', 'L')) %>%
  mutate(Group = factor(Group, c('B', 'L'))) %>%
  mutate(Group = recode(Group, B = 'A-B', L = 'A-L')) %>%
  filter(Clusters == 'CD4 T EFF') %>%
  ggplot(aes(Group, freq)) +
  geom_boxplot(fill = metacluster_colours[names(metacluster_colours) == 'CD4 T EFF'], outlier.alpha = 0, alpha = 0.5) +
  geom_jitter(position = position_jitter(0.1)) +
  theme_pubclean() +
  facet_grid(~Clusters) +
  ylab('% live cells') +
  geom_bracket(xmin = "A-B", xmax = "A-L", y.position = 10.5, inherit.aes = F,
               label = "p = 0.0034 (FDR-adjusted p < 0.01)", size = 0.5, label.size = 4, 
               tip.length = c(0.02, 0.02)) +
  theme(legend.position = 'none', 
        text = element_text(size=20),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA), 
        strip.background = element_rect(fill = "grey90", colour = "black", size = 0.5)) +
  scale_y_continuous(breaks = seq(0,20,2.5), expand = c(0.01,1))

p1_tmp

ggsave('Figures/Supplementary Plots/Unsupervised Clustering CD4 T EFF BvL v2.png', p1_tmp, width = 6, height = 6, dpi = 500)





