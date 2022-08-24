###############################################################################################
##' Manual gating MMI CyTOF data analysis
##' 21-06-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------
setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

mmi_df <- read.csv('Datatables/Manual Gating/manual_gating_mfi.csv')
rownames(mmi_df) <- NULL

new_colnames <- sub('B.cells.Naive...Median...*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'Naive B cell \\1', colnames(mmi_df))
new_colnames <- sub('B.cells.Memory...Median...*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'Memory B cell \\1', new_colnames)
new_colnames <- sub('.*NK.cells\\.\\.\\.\\.Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'NK \\1', new_colnames)
new_colnames <- sub('.*NK.cells\\.\\.\\.Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD57 +ve NK \\1', new_colnames)
new_colnames <- sub('.*CD57\\.\\.ve.Monocytes...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD57 +ve Monocyte \\1', new_colnames)
new_colnames <- sub('.*Treg...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'Treg \\1', new_colnames)
new_colnames <- sub('.*gd.T.cells...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'gd T cell \\1', new_colnames)
new_colnames <- sub('.*CD8.SEN...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD8 T SEN \\1', new_colnames)
new_colnames <- sub('.*CD8.Naive...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD8 T N \\1', new_colnames)
new_colnames <- sub('.*CD8.EMRA...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD8 T EFF \\1', new_colnames)
new_colnames <- sub('.*CD8.CM...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD8 T CM \\1', new_colnames)
new_colnames <- sub('.*CD4.SEN...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD4 T SEN \\1', new_colnames)
new_colnames <- sub('.*CD4.Naive...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD4 T N \\1', new_colnames)
new_colnames <- sub('.*CD4.EMRA...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD4 T EFF \\1', new_colnames)
new_colnames <- sub('.*CD4.CM...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'CD4 T CM \\1', new_colnames)
new_colnames <- sub('.*Non.classical.Monocytes...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'Non-classical Monocyte \\1', new_colnames)
new_colnames <- sub('.*Classical.Monocytes...Median.*(CCR2|CCR4|CCR6|CD40|CD86|Fas|PD1|CD169|CD61|CD24).*', 'Classical Monocyte \\1', new_colnames)
colnames(mmi_df) <- new_colnames

mmi_df[mmi_df == 'n/a'] <- 0

mmi_df <- mmi_df %>%
  as.data.frame() %>%
  column_to_rownames('X') %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(as.numeric)

mmi_df <- mmi_df[!rowMeans(mmi_df) < 10,]

# Differences between groups (CTRL vs ALS)
DS_res_group <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_group <- list()

for (i in seq_len(nrow(mmi_df))) {
  
  tmp_df <- data.frame(y = as.numeric(t(mmi_df[i,match(tmp_md$sample_id, colnames(mmi_df))])), tmp_md)
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_group[[i]] <- lmer(y ~ Group + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_group[[i]])
  DS_res_group[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_group) <- rownames(mmi_df)

DS_res_group <- DS_res_group %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

head(DS_res_group, 10)

# Differences between progression (F v S)
DS_res_progression <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_progression <- list()

for (i in seq_len(nrow(mmi_df))) {
  
  tmp_df <- data.frame(y = as.numeric(t(mmi_df[i,match(tmp_md$sample_id, colnames(mmi_df))])), tmp_md)
  tmp_df <- tmp_df %>%
    dplyr::filter(Group == 'ALS')
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_progression[[i]] <- lmer(y ~ Progression + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_progression[[i]])
  DS_res_progression[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_progression) <- rownames(mmi_df)

DS_res_progression <- DS_res_progression %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

head(DS_res_progression,10)

# Differences between onset (F v S)
DS_res_onset <- data.frame('Chisq' = NA, 'Df' = NA, 'p_val' = NA)
lm1_ds_onset <- list()

for (i in seq_len(nrow(mmi_df))) {
  
  tmp_df <- data.frame(y = as.numeric(t(mmi_df[i,match(tmp_md$sample_id, colnames(mmi_df))])), tmp_md)
  tmp_df <- tmp_df %>%
    dplyr::filter(Group == 'ALS')
  
  lm0 <- lmer(y ~ 1 + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  lm1_ds_onset[[i]] <- lmer(y ~ Onset + (1|CyTOF_run_date), data = tmp_df, weights = tmp_df$col_sums)
  
  res <- anova(lm0, lm1_ds_onset[[i]])
  DS_res_onset[i,] <- cbind(as.numeric(c(res$Chisq[2], res$Df[2], res$`Pr(>Chisq)`[2])))
}

rownames(DS_res_onset) <- rownames(mmi_df)

DS_res_onset <- DS_res_onset %>%
  mutate_all(function(x) round(x, 4)) %>%
  rownames_to_column('Metacluster') %>%
  arrange(p_val)

head(DS_res_onset, 10)

# Write results to folder
write.csv(DS_res_group, 'Datatables/Results/Manual Gating Analysis/Differential States/LMER Group CTRLvALS results.csv', row.names = F)
write.csv(DS_res_progression, 'Datatables/Results/Manual Gating Analysis/Differential States/LMER Progression SvF results.csv', row.names = F)
write.csv(DS_res_onset, 'Datatables/Results/Manual Gating Analysis/Differential States/LMER Onset BvL results.csv', row.names = F)

