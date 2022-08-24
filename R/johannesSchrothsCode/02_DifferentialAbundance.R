###############################################################################################
##' CyTOF data analysis of ALS patient PBMCs
##' 29-03-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------------------------------------------------

#rm(list = ls())

set.seed(1234)

setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')
library(dplyr)
library(edgeR)

## Differential abundance ---------------------------------------------------------------------------------------------------------------

counts <- df %>%
  dplyr::count(metaclusters, sample_id) %>%
  spread(key = 'sample_id', value = 'n') %>%
  filter(metaclusters != 'Excluded') %>%
  column_to_rownames('metaclusters') %>%
  as.data.frame.matrix() %>%
  replace(is.na(.), 0)

counts <- counts[!rowMeans(counts) < 5,]

group <- as.factor(md[match(colnames(counts), md$sample_id), 'Group'])
progression <- as.factor(md[match(colnames(counts), md$sample_id), 'Progression'])
onset <- as.factor(md[match(colnames(counts), md$sample_id), 'Onset'])
riluzole <- as.factor(md[match(colnames(counts), md$sample_id), 'Riluzole'])
date <- as.factor(md[match(colnames(counts), md$sample_id), 'CyTOF_run_date'])
nrow(counts)
DA_res_group <- counts %>% 
  DGEList(lib.size=colSums(counts)) %>%
  estimateDisp(design = model.matrix(~ 0 + group + date)) %>%
  glmQLFit(robust = T) %>%
  glmQLFTest(fit, contrast = makeContrasts(CTRLvsALS = groupCTRL - groupALS, 
                                           levels = model.matrix(~ 0 + group + date))) %>%
  .$table %>%
  mutate(FDR = p.adjust(PValue, 'fdr')) %>%
  rownames_to_column('Celltype') %>%
  arrange(FDR)

DA_res_progression <- counts %>% 
  DGEList(lib.size=colSums(counts)) %>%
  estimateDisp(design = model.matrix(~ 0 + progression + date)) %>%
  glmQLFit(robust=F) %>%
  glmQLFTest(fit, contrast = makeContrasts(FvsS = progressionF - progressionS, 
                                           levels = model.matrix(~ 0 + progression + date))) %>%
  .$table %>%
  mutate(FDR = p.adjust(PValue, 'fdr')) %>%
  rownames_to_column('Celltype') %>%
  arrange(FDR)

DA_res_onset <- counts %>% 
  DGEList(lib.size=colSums(counts)) %>%
  estimateDisp(design = model.matrix(~ 0 + onset + date)) %>%
  glmQLFit(robust=TRUE) %>%
  glmQLFTest(fit, contrast = makeContrasts(BvsL = onsetB - onsetL, 
                                           levels = model.matrix(~ 0 + onset + date))) %>%
  .$table %>%
  mutate(FDR = p.adjust(PValue, 'fdr')) %>%
  rownames_to_column('Celltype') %>%
  arrange(FDR)


DA_res_riluzole <- counts %>% 
  DGEList(lib.size=colSums(counts)) %>%
  estimateDisp(design = model.matrix(~ 0 + riluzole)) %>%
  glmQLFit(robust=TRUE) %>%
  glmQLFTest(fit, contrast = makeContrasts(YesvNo = riluzoleyes - riluzoleno,
                                           levels = model.matrix(~ 0 + riluzole))) %>%
  .$table %>%
  mutate(FDR = p.adjust(PValue, 'fdr')) %>%
  rownames_to_column('Celltype') %>%
  arrange(FDR)



# Write results to folder
write.csv(DA_res_group, 'Datatables/Results/Differential Abundance/Group CTRLvALS results.csv', row.names = F)
write.csv(DA_res_progression, 'Datatables/Results/Differential Abundance/Progression SvF results.csv', row.names = F)
write.csv(DA_res_onset, 'Datatables/Results/Differential Abundance/Onset BvL results.csv', row.names = F)
write.csv(DA_res_riluzole, 'Datatables/Results/Differential Abundance/Riluzole YvN results.csv', row.names = F)

