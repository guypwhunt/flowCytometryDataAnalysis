###############################################################################################
##' Manual gating percent positive CyTOF data analysis
##' 21-06-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------
setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')


## Clinical correlation with cluster abundance -----------------------------------------------------------------------------------------

variables <- rownames(counts)[-grep('Undefined', rownames(counts))]

clinical_parameters <- cbind(md[match(colnames(counts), md$sample_id),], t(counts)) %>%
  dplyr::filter(Group == 'ALS') %>%
  dplyr::filter(!sample_id == 'BLT00092') %>%
  dplyr::select(-grep('date', colnames(md)), -sample_id, -file_name, -Group, -Ethnicity, -DOB,
                -Lactate, -CK, -NfH_at_V1, -Disease_duration_to_death_months, -C9orf72, -batch, -type, -Progression) %>%
  mutate_if(is.factor, as.integer) %>%
  mutate_all(as.numeric)

a <- clinical_parameters[,colnames(clinical_parameters) %in% colnames(md)]
b <- clinical_parameters[,colnames(clinical_parameters) %in% rownames(counts)]
c <- md[match(colnames(counts), md$sample_id) & md$Group == 'ALS' & md$sample_id != 'BLT00092', 'CyTOF_run_date']
d <- md[match(colnames(counts), md$sample_id) & md$Group == 'ALS' & md$sample_id != 'BLT00092', 'sample_id']
w <- tmp_md[match(colnames(counts), tmp_md$sample_id) & tmp_md$Group == 'ALS' & tmp_md$sample_id != 'BLT00092', 'col_sums']

res_clin_corr <- data.frame(matrix(ncol = ncol(b), nrow = ncol(a)))
lm1_clin_corr <- list()
clin_corr_names <- vector()


for (i in 1:ncol(a)) {
  for(j in 1:ncol(b)) {
    tmp_df <- data.frame('y' = b[,j], 'x' = a[,i], 'batch' = c)
    
    lm0 <- lme4::lmer(y ~ 1 + (1|c), REML = F,  tmp_df)
    lm1_clin_corr[[paste(colnames(b)[j], colnames(a)[i])]] <- lme4::lmer(y ~ x + (1|c),  REML = F, tmp_df)
    res <- anova(lm0, lm1_clin_corr[[paste(colnames(b)[j], colnames(a)[i])]])
    res_clin_corr[i,j] <- res$`Pr(>Chisq)`[2]
  }
}

rownames(res_clin_corr) <- colnames(a)
colnames(res_clin_corr) <- colnames(b)

clin_corr_signif <- res_clin_corr %>%
  as.data.frame() %>%
  rownames_to_column('clinical_parameter') %>%
  gather(key = 'Metacluster', value = 'p_val', -clinical_parameter) %>%
  dplyr::filter(p_val < 0.05)

clin_corr_signif %>%
  arrange(p_val)

clin_par_signif <- clinical_parameters %>%
  as.data.frame() %>%
  dplyr::select(clin_corr_signif$clinical_parameter, clin_corr_signif$Metacluster)

effects_df <- as.data.frame(effects::effect(term = 'x', mod = lm1_clin_corr$`CD8 T SEN ALSFRSR_at_V1`))

performance::check_model(lm1_clin_corr$`CD8 T SEN ALSFRSR_at_V1`)

res2 <-cor.test(clin_par_signif$Age_at_V1, clin_par_signif$`CD8 T EFF`,  method = "spearman")
res2

shapiro.test(clin_par_signif$`CD8 T CM`)

corrr <- cor(na.omit(clinical_parameters), method = 'spearman')


p.mat <- cor_pmat(na.omit(clinical_parameters), method = 'spearman')
ggcorrplot(corrr, p.mat = p.mat, insig = 'blank', type = 'lower',sig.level = 0.05)

clin_par_signif %>%
  ggplot(aes(`CD8 T CM`, PRL)) +
  geom_point() +
  ggpubr::stat_cor(method = 'spearman')


ggplot() +
  geom_point(data = clin_par_signif, aes(ALSF, `Naive B cell`)) +
  geom_line(data=effects_df, aes(x=x, y=fit), color="blue") +
  geom_ribbon(data=effects_df, aes(x=x, ymin=lower, ymax=upper), alpha= 0.3, fill="blue")

lm1_clin_corr$`CD8 T SEN Months_onset_to_V1`


MuMIn::r.squaredGLMM(lm1_clin_corr$`CD8 T EFF Age_at_V1`)


clin_corr_signif

clin_corr_plots <- list()

for(i in 1:nrow(clin_corr_signif)) {
  
  tmp_df <- clin_par_signif %>%
    dplyr::select(clin_corr_signif$clinical_parameter[i], clin_corr_signif$Metacluster[i])
  
  colnames(tmp_df) <- c('x', 'y')
  
  clin_corr_plots[[i]] <- tmp_df %>%
    ggplot(aes(x,y)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F) +
    ggpubr::stat_cor(method = 'kendall') +
    theme_pubclean() +
    labs(x = clin_corr_signif$clinical_parameter[i],
         y = clin_corr_signif$Metacluster[i]) +
    xlim(0,NA) +
    ylim(0,NA)
  
}

clin_corr_plots[[4]]
