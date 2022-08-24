###############################################################################################
##' Manual gating percent positive CyTOF data analysis
##' 21-06-2021 Johannes Schroth
###############################################################################################

## Setup --------------------------------------------------------------------------------------
setwd('~/Dropbox/PhD/01 Data/Other/Ozlem/CyTOF data/Johannes CyTOF Analysis/Combined/')

### Survival Analysis -------------------------------------------------------------------------------------------------------------------

counts_batch_corrected <- limma::removeBatchEffect(counts, md[match(colnames(counts), md$sample_id),'CyTOF_run_date'],
                                                   design = model.matrix(~md[match(colnames(counts), md$sample_id),'Group']))


surv_data <- cbind(md[match(colnames(counts_batch_corrected), md$sample_id),], t(counts_batch_corrected)) %>%
  as.data.frame() %>%
  dplyr::filter(Group == 'ALS') %>%
  tibble(.name_repair = 'universal')

variables <- t(counts_batch_corrected) %>%
  as.data.frame() %>%
  tibble(.name_repair = 'universal') %>%
  colnames

res_cox <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(res_cox) <- c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
cox_plots <- c()

surv_data_tmp <- surv_data %>%
  dplyr::select(variables) %>%
  mutate_all(funs(cut(., breaks = c(-Inf, median(.), Inf), labels = c('low', 'high')))) %>%
  cbind('Survival' = surv_data$Survival_from_onset_months, 
        'Outcome' = surv_data$Outcome)  %>%
  tibble(.name_repair = 'universal')

for (i in seq(variables)) {
  
  univ_model <- coxph(as.formula(paste0('Surv(Survival, Outcome) ~', '`', variables[i], '`')), data = surv_data_tmp)
  
  x <- summary(univ_model)
  res_cox[variables[i],] <- c(signif(x$coef[1], digits=3),
                              paste0(signif(x$coef[2], digits=2), " (", 
                                     signif(x$conf.int[,"lower .95"], 2), "-", 
                                     signif(x$conf.int[,"upper .95"],2), ")"),
                              signif(x$wald["test"], digits=3),
                              signif(x$wald["pvalue"], digits=3))
  
}

res_cox <- res_cox %>%
  rownames_to_column('Metacluster') %>%
  mutate(Metacluster = gsub('.ve', '+ve', Metacluster)) %>%
  mutate(Metacluster = gsub('\\.', ' ', Metacluster)) %>%
  arrange(p.value)
res_cox

write.csv(res_cox, 'Datatables/Results/Cox Survival Model Results.csv')

ggsurvplot(survfit(Surv(surv_data_tmp$Survival, surv_data_tmp$Outcome) ~ CD4.T.SEN, data = surv_data_tmp),
           pval=T, conf.int=T) +
  ggtitle('CD4 T senescent cells')
