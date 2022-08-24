###############################################################################################
##' Manual gating percent positive CyTOF data analysis
##' 21-06-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------
setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

count_df <- read.csv('Datatables/Manual Gating/manual_gating_count_df.csv')

count_df <- count_df %>%
  as.data.frame() %>%
  column_to_rownames('X')

colnames(count_df) <- c('Memory B cell', 'Naive B cell', 'CD57 +ve NK cells', 'NK Cells', 'CD57 +ve Monocytes',
                        'Classical Monocytes', 'Non-classical Monocytes', 'Treg', 'CD4 T SEN', 'CD4 T CM', 'CD4 T EFF',
                        'CD4 T N', 'CD8 T SEN', 'CD8 T CM', 'CD8 T EFF', 'CD8 T N', 'gd T cell')

count_df <- count_df %>%
  t()


## Analysis ------------------------------------------------------------------------------------
count_df <- count_df[!rowMeans(count_df) < 5,]

group <- as.factor(md[match(colnames(count_df), md$sample_id), 'Group'])
progression <- as.factor(md[match(colnames(count_df), md$sample_id), 'Progression'])
onset <- as.factor(md[match(colnames(count_df), md$sample_id), 'Onset'])
riluzole <- as.factor(md[match(colnames(count_df), md$sample_id), 'Riluzole'])
date <- as.factor(md[match(colnames(count_df), md$sample_id), 'CyTOF_run_date'])

DA_res_group <- count_df %>% 
  DGEList(lib.size=colSums(count_df)) %>%
  estimateDisp(design = model.matrix(~ 0 + group + date)) %>%
  glmQLFit(robust = T) %>%
  glmQLFTest(fit, contrast = makeContrasts(CTRLvsALS = groupCTRL - groupALS, 
                                           levels = model.matrix(~ 0 + group + date))) %>%
  .$table %>%
  rownames_to_column('Celltype') %>%
  arrange(PValue)

DA_res_progression <- count_df %>% 
  DGEList(lib.size=colSums(count_df)) %>%
  estimateDisp(design = model.matrix(~ 0 + progression + date)) %>%
  glmQLFit(robust=F) %>%
  glmQLFTest(fit, contrast = makeContrasts(FvsS = progressionF - progressionS, 
                                           levels = model.matrix(~ 0 + progression + date))) %>%
  .$table %>%
  rownames_to_column('Celltype') %>%
  arrange(PValue)

DA_res_onset <- count_df %>% 
  DGEList(lib.size=colSums(count_df)) %>%
  estimateDisp(design = model.matrix(~ 0 + onset + date)) %>%
  glmQLFit(robust=TRUE) %>%
  glmQLFTest(fit, contrast = makeContrasts(BvsL = onsetB - onsetL, 
                                           levels = model.matrix(~ 0 + onset + date))) %>%
  .$table %>%
  rownames_to_column('Celltype') %>%
  arrange(PValue)

DA_res_riluzole <- count_df %>% 
  DGEList(lib.size=colSums(count_df)) %>%
  estimateDisp(design = model.matrix(~ 0 + riluzole)) %>%
  glmQLFit(robust=TRUE) %>%
  glmQLFTest(fit, contrast = makeContrasts(YesvNo = riluzoleyes - riluzoleno,
                                           levels = model.matrix(~ 0 + riluzole))) %>%
  .$table %>%
  rownames_to_column('Celltype') %>%
  arrange(PValue)

# Write results to folder
write.csv(DA_res_group, 'Datatables/Results/Manual Gating Analysis/Differential Abundance/Group CTRLvALS results manual gating.csv', row.names = F)
write.csv(DA_res_progression, 'Datatables/Results/Manual Gating Analysis/Differential Abundance/Progression SvF results manual gating.csv', row.names = F)
write.csv(DA_res_onset, 'Datatables/Results/Manual Gating Analysis/Differential Abundance/Onset BvL results manual gating.csv', row.names = F)
write.csv(DA_res_riluzole, 'Datatables/Results/Manual Gating Analysis/Differential Abundance/Riluzole YvN results manual gating.csv', row.names = F)
