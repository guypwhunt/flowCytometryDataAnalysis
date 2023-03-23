library(vctrs)
library(limma)
library(umap)
library(readxl)
library(data.table)
library(dplyr)
library(ggrepel)
library(tibble)
library(factoextra)
library(tidyr)
library(ltm)

dir.create("data/lipidomics/Results/pcaAnalysis", showWarnings = FALSE)

visit <- 12

sampleInformation <-
  read_excel("data/metadata/clinicalData.xlsx") %>%
  as.data.frame() %>%
  #filter(visit == visit | caseControl != "Case")  %>%
  filter(experiment == "lipidomics") %>%
  dplyr::select(classification, patient_id, sample_id, visit, caseControl, earlyLate,
         caseControlExperiment, sightOnsetExperiment, progressionExperiment,
         progressionAndSightOnsetExperiment, fastSlow, fastSlow_original, BulbarLimb,
         gender, ethnicity, statinUse, timeFromVisit1InYears,
         ageAtVisit, ageAtOnset, diagnosticDelayInYears, diseaseDurationInYears,
         sampleStorageDays, timeFromOnsetToVisitInYears, alsfrsR,
         survivalFromVisit, delataAlsfrsRScore
         )

expressionData <- fread(
  #"data/lipidomics/normalisedasinhTransformedExpressionDataRawOutliersAndDuplicatesRemoved.csv"
  "data/lipidomics/expressionDataRaw.csv"
  ) %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  asinh() %>%
  normalizeBetweenArrays() #%>%
  # dplyr::filter(grepl("Rv",V1)) %>%
  # dplyr::filter(nchar(V1)==4) %>%


sampleInformation <- sampleInformation %>%
  filter(classification %in% colnames(expressionData))

expressionData <- expressionData[, unique(sampleInformation$classification)]

pcaResults <- expressionData %>%
  na.omit() %>%
  prcomp()

eigenValue <- get_eigenvalue(pcaResults)

variableStats <- get_pca_var(pcaResults)

rotation <-
  unclass(pcaResults$rotation) %>% as.data.frame() %>%
  rownames_to_column("ID") %>%
  merge(sampleInformation, by.x = "ID", by.y = "classification")

pcaIndividualsResults <- as.data.frame(pcaResults$x) %>%
  rownames_to_column("ID")

pcaIndividualsResults$group <- NA

spmNames <- unique(pcaIndividualsResults$ID)

#Add Group Information
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("Rv", spmNames)], "group"] <- "Resolvins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("PD", spmNames)], "group"] <- "Protectins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("diHDHA", spmNames)], "group"] <- "Maresins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("MaR1", spmNames)], "group"] <- "Maresins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("LX", spmNames)], "group"] <- "Lipoxins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("5S,15,S-diHETE", spmNames)], "group"] <- "Lipoxins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("LT", spmNames)], "group"] <- "Leukotriene B4 Metabolome"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("5S,12S-diHETE", spmNames)], "group"] <- "Leukotriene B4 Metabolome"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("PG", spmNames)], "group"] <- "Prostaglandins"
pcaIndividualsResults[pcaIndividualsResults$ID %in% spmNames[grep("Tx", spmNames)], "group"] <- "Thromboxane B2"

# Individuals Plot
ggplot(pcaIndividualsResults, aes(x=PC1, y=PC2, colour = group, label = ID)) +
  geom_point() +
  geom_label_repel()

rotation$statinUse <- factor(rotation$statinUse, levels = c("Yes", "No"))

# Variables Plot
ggplot(rotation, aes(x=PC1, y=PC2, colour = caseControl)) +
  geom_point() +
  guides(color = guide_legend(title = "Statin Use"))

# ## Experiment Variables Plot
# for (column in colnames(rotation)){
#   try({
#     print({
#       ggplot(rotation, aes_string(x = "PC3", y = "PC4", colour = column)) +
#         geom_point()
#     })
#   }
#   )
# }

continousDFCorrelations <- rotation %>%
  select_if(is.numeric) %>%
  as.matrix() %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value,-var1) %>%
  arrange(value) %>%
  na.omit()

categoricalDF <- rotation %>%
  #select_if(is.character) %>%
  dplyr::select(caseControl, earlyLate, fastSlow, BulbarLimb,
         gender, statinUse) %>%
  cbind(rotation[, grepl("PC", colnames(rotation))])

categoricalDFColNames <- colnames(categoricalDF)

categoricalDFCorrelations <- data.frame(var1 = as.character(),
                                        var2 = as.character(),
                                        value = as.numeric())

for(clinicalColumnName in categoricalDFColNames[!grepl("PC", categoricalDFColNames)]) {
  for (pcColumnName in categoricalDFColNames[grepl("PC", categoricalDFColNames)]) {
    tempDf <- categoricalDF[, c(clinicalColumnName, pcColumnName)] %>%
      na.omit()

    try({
      categoricalDFCorrelations[nrow(categoricalDFCorrelations) + 1,] <-
        c(clinicalColumnName,
          pcColumnName,
          biserial.cor(tempDf[, pcColumnName],
                       tempDf[, clinicalColumnName]))
    })
  }
}

correlationResults <- categoricalDFCorrelations %>%
  rbind(continousDFCorrelations) %>%
  na.omit()

for(variable in c("var1", "var2")) {
  correlationResults[correlationResults[, variable] == 'ageAtOnset' ,variable] <- 'Age at Onset (Years)'
  correlationResults[correlationResults[, variable] == 'ageAtVisit' ,variable] <- 'Age at Visit (Years)'
  correlationResults[correlationResults[, variable] == 'alsfrsR' ,variable] <- 'ALSFRS'
  correlationResults[correlationResults[, variable] == 'caseControl' ,variable] <- 'Case Control'
  correlationResults[correlationResults[, variable] == 'delataAlsfrsRScore' ,variable] <- 'Delta ALSFRS'
  correlationResults[correlationResults[, variable] == 'diagnosticDelayInYears' ,variable] <- 'Diagnostic Delay (Years)'
  correlationResults[correlationResults[, variable] == 'diseaseDurationInYears' ,variable] <- 'Disease Duration (Years)'
  correlationResults[correlationResults[, variable] == 'fastSlow' ,variable] <- 'Progression'
  correlationResults[correlationResults[, variable] == 'sampleStorageDays' ,variable] <- 'Sample Storage Time (Days)'
  correlationResults[correlationResults[, variable] == 'gender' ,variable] <- 'Sex'
  correlationResults[correlationResults[, variable] == 'BulbarLimb' ,variable] <- 'Site of Onset'
  correlationResults[correlationResults[, variable] == 'statinUse' ,variable] <- 'Statin'
  correlationResults[correlationResults[, variable] == 'survivalFromVisit' ,variable] <- 'Survival From Visit (Years)'
  correlationResults[correlationResults[, variable] == 'timeFromVisit1InYears' ,variable] <- 'Time from First Visit (Years)'
  correlationResults[correlationResults[, variable] == 'timeFromOnsetToVisitInYears' ,variable] <- 'Time from Onset to Visit (Years)'
  correlationResults[correlationResults[, variable] == 'earlyLate' ,variable] <- 'Visit'
  correlationResults[correlationResults[, variable] == 'visit' ,variable] <- 'Visit Number'
}

colnames(correlationResults) <- c("Variable 1", "Variable 2","Correlation")

fwrite(correlationResults, paste0("data/lipidomics/Results/pcaAnalysis/pcaCorrelationResultsVisit", visit, ".csv"))
