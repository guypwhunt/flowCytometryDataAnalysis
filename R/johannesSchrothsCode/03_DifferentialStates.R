###############################################################################################
##' CyTOF data analysis of ALS patient PBMCs
##' 29-03-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------

#rm(list = ls())

set.seed(1234)

setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

## Differential states ------------------------------------------------------------------------

# Calculate median expression of functional markers for each metacluster celltype
medians <- df %>%
  dplyr::filter(!metaclusters %in% c('Excluded', 'Undefined')) %>%
  group_by(sample_id, metaclusters) %>%
  dplyr::select(sample_id, metaclusters, panel$antigen[panel$marker_class == 'state']) %>%
  summarise_if(is.numeric, median) %>%
  gather(key = 'key', value = 'value', -sample_id, -metaclusters) %>%
  spread(key = 'sample_id', value = 'value') %>%
  unite('rows', metaclusters:key) %>%
  column_to_rownames('rows') %>%
  replace(., is.na(.), 0)


# Remove cells with no/low expression of marker
medians <- medians[!rowMeans(medians) < 0.1,]

# Differences between groups (CTRL vs ALS)
DS_res_group <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_group <- list()

for (i in seq_len(nrow(medians))) {
  
  tmp_df <- data.frame(y = as.numeric(t(medians[i,match(tmp_md$sample_id, colnames(medians))])), tmp_md)
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_group[[i]] <- lmer(y ~ Group + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_group[[i]])
  DS_res_group[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_group) <- rownames(medians)

DS_res_group <- DS_res_group %>%
  mutate(p_adj = p.adjust(p_val, method = 'fdr')) %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

DS_res_group

# Differences between progression (F v S)
DS_res_progression <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_progression <- list()

for (i in seq_len(nrow(medians))) {
  
  tmp_df <- data.frame(y = as.numeric(t(medians[i,match(tmp_md$sample_id, colnames(medians))])), tmp_md)
  tmp_df <- tmp_df %>%
    dplyr::filter(Group == 'ALS')
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_progression[[i]] <- lmer(y ~ Progression + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_progression[[i]])
  DS_res_progression[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_progression) <- rownames(medians)

DS_res_progression <- DS_res_progression %>%
  mutate(p_adj = p.adjust(p_val, method = 'fdr')) %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

DS_res_progression

# Differences between onset (F v S)
DS_res_onset <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_onset <- list()

for (i in seq_len(nrow(medians))) {
  
  tmp_df <- data.frame(y = as.numeric(t(medians[i,match(tmp_md$sample_id, colnames(medians))])), tmp_md)
  tmp_df <- tmp_df %>%
    dplyr::filter(Group == 'ALS')
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_onset[[i]] <- lmer(y ~ Onset + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_onset[[i]])
  DS_res_onset[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_onset) <- rownames(medians)

DS_res_onset <- DS_res_onset %>%
  mutate(p_adj = p.adjust(p_val, method = 'fdr')) %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

head(DS_res_onset, 10)

# Write results to folder
write.csv(DS_res_group, 'Datatables/Results/Differential States/LMER Group CTRLvALS results.csv', row.names = F)
write.csv(DS_res_progression, 'Datatables/Results/Differential States/LMER Progression SvF results.csv', row.names = F)
write.csv(DS_res_onset, 'Datatables/Results/Differential States/LMER Onset BvL results.csv', row.names = F)

