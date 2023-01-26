library(limma)
library(umap)
library(dplyr)
library(ggplot2)

figureDirectory <- "data/lipidomics/Results/"

fileNames <- list.files(figureDirectory, pattern =".csv")

fileNamesPath <- paste0(figureDirectory, fileNames)

dfs <- lapply(fileNamesPath, read.csv)

names(dfs) <- fileNames

i <- 1

for(df in dfs){
  df$experiment <- names(dfs)[i]
  i <- i + 1

  if (exists("combinedDf")){
    combinedDf <- rbind(combinedDf, df)
  } else {
    combinedDf <- df
  }
}

unique(combinedDf$experiment)

combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowFastLateVsFastEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupFastLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastLimbLate.csv",]
#combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "BulbarVsLimbBulbarLateVsBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowSlowEarlyVsFastEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowSlowLateVsFastLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "EarlyVsLateLateVsEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowFastLateVsSlowEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowSlowEarlyVsFastLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupSlowLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "BulbarVsLimbBulbarLateVsLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "BulbarVsLimbLimbEarlyVsBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupFastBulbarLate.csv",]
combinedDf <- combinedDf[!combinedDf$experiment %in%
                           c("BulbarVsLimbBulbarEarlyVsLimbLate.csv",
                             "BulbarVsLimbLimbEarlyVsBulbarEarly.csv",
                             "BulbarVsLimbLimbLateVsBulbarEarly.csv",
                             "BulbarVsLimbLimbLateVsBulbarLate.csv",
                             "BulbarVsLimbLimbLateVsLimbEarly.csv",
                             "FastVsSlowFastEarlyVsSlowLate.csv",
                             "FastVsSlowSlowEarlyVsSlowLate.csv",
                             "FastVsSlowSlowLateVsFastEarly.csv",
                             "FastVsSlowSlowLateVsSlowEarly.csv"
                           ),]


combinedDf <- combinedDf[!grepl("ControlVs", combinedDf$experiment, fixed = TRUE),]



combinedDf[combinedDf$experiment == "FastVsSlowFastEarlyVsFastLate.csv", "experiment"] <- "Early Fast vs Late Fast"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowBulbarLate.csv", "experiment"] <- "Late Fast Bulbar vs Late Slow Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastLimbLate.csv", "experiment"] <- "Late Slow Limb vs Late Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowLimbLate.csv","experiment"] <- "Late Fast Limb vs Late Slow Limb"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbLateVsControl.csv","experiment"] <- "Late Limb vs Controls"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbEarlyVsControl.csv","experiment"] <- "Early Limb vs Controls"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlBulbarLateVsControl.csv","experiment"] <- "Late Bulbar vs Controls"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlBulbarEarlyVsControl.csv","experiment"] <- "Early Bulbar vs Controls"
combinedDf[combinedDf$experiment == "ProgressionVsControlSlowLateVsControl.csv","experiment"] <- "Late Slow vs Controls"
combinedDf[combinedDf$experiment == "ProgressionVsControlSlowEarlyVsControl.csv","experiment"] <- "Early Slow vs Controls"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastLateVsControl.csv","experiment"] <- "Late Fast vs Controls"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastEarlyVsControl.csv","experiment"] <- "Early Fast vs Controls"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastEarlyVsControl.csv","experiment"] <- "Early Fast vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlFastLimbLateVsControl.csv","experiment"] <- "Late Fast Limb vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowBulbarEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Slow Bulbar"
#combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastBulbarEarly.csv","experiment"] <- "Early Fast Limb vs Early Fast Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastLimbLate.csv","experiment"] <- "Early Fast Limb vs Late Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowBulbarEarlyVsControl.csv","experiment"] <- "Early Slow Bulbar vs Controls"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarEarlyVsBulbarLate.csv","experiment"] <- "Early Bulbar vs Late Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupSlowBulbarLate.csv","experiment"] <- "Early Slow Bulbar vs Late Slow Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowLimbLate.csv" ,"experiment"] <- "Late Slow Bulbar vs Late Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowLimbLate.csv" ,"experiment"] <- "Late Slow Bulbar vs Late Slow Limb"
combinedDf[combinedDf$experiment == "FastVsSlowFastEarlyVsSlowEarly.csv" ,"experiment"] <- "Early Fast vs Early Slow"
combinedDf[combinedDf$experiment == "FastVsSlowFastLateVsSlowLate.csv" ,"experiment"] <- "Late Fast vs Late Slow"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupSlowLimbEarly.csv" ,"experiment"] <- "Early Slow Bulbar vs Early Slow Limb"
combinedDf[combinedDf$experiment == "CaseVsControlEarlyVsControl.csv" ,"experiment"] <- "Early vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastBulbarLate.csv" ,"experiment"] <- "Early Fast Bulbar vs Late Fast Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupFastLimbLate.csv" ,"experiment"] <- "Late Fast Bulbar vs Late Fast Limb"
combinedDf[combinedDf$experiment == "EarlyVsLateEarlyVsLate.csv" ,"experiment"] <- "Early vs Late"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupSlowLimbLate.csv" ,"experiment"] <- "Early Slow Limb vs Late Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowLimbEarly.csv" ,"experiment"] <- "Early Fast Limb vs Early Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowBulbarLateVsControl.csv" ,"experiment"] <- "Slow Bulbar vs Controls"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarEarlyVsLimbEarly.csv" ,"experiment"] <- "Early Bulbar vs Early Limb"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarLateVsLimbLate.csv" ,"experiment"] <- "Late Bulbar vs Late Limb"
combinedDf[combinedDf$experiment == "BulbarVsLimbLimbEarlyVsLimbLate.csv" ,"experiment"] <- "Early Limb vs Late Limb"
combinedDf[combinedDf$experiment == "CaseVsControlLateVsControl.csv" ,"experiment"] <- "Late vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlFastBulbarEarlyVsControl.csv" ,"experiment"] <- "Early Fast Bulbar vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlFastBulbarLateVsControl.csv" ,"experiment"] <- "Late Fast Bulbar vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlFastLimbEarlyVsControl.csv" ,"experiment"] <- "Early Fast Limb vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowLimbEarlyVsControl.csv" ,"experiment"] <- "Early Slow Limb vs Controls"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowLimbLateVsControl.csv" ,"experiment"] <- "Late Slow Limb vs Controls"


colnames(combinedDf)

combinedDf <- combinedDf[, c(
  "X",
  "experiment",
  "logFC",
  "AveExpr",
  "t",
  "B",
  "P.Value",
  "adj.P.Val"
)]

combinedDf$logPValue = 0-log10(combinedDf$P.Value)
combinedDf$logFDRPValue = 0-log10(combinedDf$adj.P.Val)


colnames(combinedDf) <- c(
  "Lipid",
  "Comparison",
  "Log Fold Change",
  "Average Expression",
  "T Statistic",
  "B Statistic (Log Odds)",
  "Raw P-Value",
  "FDR Adjusted P-Value",
  "Minus Log (Raw P-Value)",
  "Minus Log (FDR Adjusted P-Value)"
)

dir.create("data/pValueAdjustmentsResults", showWarnings = FALSE)

fwrite(combinedDf,
       paste0("data/pValueAdjustmentsResults/lipidomics.csv"))
