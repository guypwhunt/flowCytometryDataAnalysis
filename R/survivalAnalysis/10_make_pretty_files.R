library(dplyr)
library(ggplot2)

figureDirectory <- "data/survivalAnalysis/"

fileNames <- list.files(figureDirectory, pattern =".csv")

fileNames <- fileNames[grepl("DiseaseDurationBiologicalModel", fileNames) | grepl("DiseaseDurationClinicalModel", fileNames)]

fileNamesPath <- paste0(figureDirectory, fileNames)

dfs <- lapply(fileNamesPath, read.csv)

names(dfs) <- fileNames

i <- 1

variableNames <- c()

for (df in dfs) {
  colnames(df) <- c("Variable", "Regression Coefficient", "Hazard Ratio",
                    "Confidence Interval", "Statistical Significance", "Global Significance (P-Value)")


  df[df$Variable == "ethnicityID2", "Variable"] <- "Ethnicity (Caucasian)"
  df[df$Variable == "ethnicityID3", "Variable"] <- "Ethnicity (Indian)"
  df[df$Variable == "alsfrsR", "Variable"] <- "ALSFRS-R Score"
  df[df$Variable == "riluzoleyes", "Variable"] <- "Riluzole"
  df[df$Variable == "ageAtOnset", "Variable"] <- "Age at Onset (Years)"
  df[df$Variable == "delataAlsfrsRScore", "Variable"] <- "ALSFRS-R Score Progression Rate"
  df[df$Variable == "GPR32_median_Positive_Naïve_B_Cells", "Variable"] <- "N-B (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_HLA_Negative_DR_Negative__Classical_Monocytes", "Variable"] <- "HLA-DR-C-M (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Classical_Monocytes_.CD11b_Low.", "Variable"] <- "CD11b Low HLA-DR-aCD11b+C-M (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_HLA_Negative_DR_Negative__Intermediate_Monocytes", "Variable"] <- "HLA-DR-IM-M (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_HLA_Negative_DR_Negative__Activated_CD11b_Positive__Intermediate_Monocytes_.CD11b_Low.", "Variable"] <- "CD11b Low HLA-DR-aCD11b+IM-M (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_Follicular_B_Cells", "Variable"] <- "F-B (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Positive_.", "Variable"] <- "CD24+US-B (Median GPR32)"
  df[df$Variable == "GPR32_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.", "Variable"] <- "CD24-US-B (Median GPR32)"

  df[df$Variable == "GPR18_median_Positive_Follicular_B_Cells", "Variable"] <- "F-B (Median GPR18)"
  df[df$Variable == "GPR18_median_Positive_Unswitched_Memory_B_Cells_.CD24_Negative_.", "Variable"] <- "CD24-US-B (Median GPR18)"

  variableNames <- append(variableNames, df$Variable)

  fwrite(df, paste0("data/pValueAdjustmentsResults/", names(dfs[i])))
  i <<- i + 1
}

unique(variableNames)


df <- read.csv(
  paste0(
    "data/survivalAnalysis/compareAllModels/",
    "diseaseDuration.csv"
  ),
  row.names = 1
)

df$Model <- rownames(df)
colnames(df) <- c("Log-Likelihood",
                  "Chi-squared",
                  "Degrees of Freedom",
                  "Global Significance (P-Value)",
                  "Model")

df <- df[, c("Model",
             "Log-Likelihood",
             "Chi-squared",
             "Degrees of Freedom",
             "Global Significance (P-Value)")]

df[df$Model == 1, "Model"] <- "Clinical"
df[df$Model == 2, "Model"] <- "GPR18 (FlowSOM)"
df[df$Model == 3, "Model"] <- "GPR18 (Phenograph)"
df[df$Model == 4, "Model"] <- "GPR32 (FlowSOM)"
df[df$Model == 5, "Model"] <- "GPR32 (Phenograph)"
df[df$Model == 6, "Model"] <- "GPR18 and GPR32 (FlowSOM)"
df[df$Model == 7, "Model"] <- "GPR18 and GPR32 (Phenograph)"

fwrite(df, paste0("data/pValueAdjustmentsResults/", "diseaseDuration.csv"))
