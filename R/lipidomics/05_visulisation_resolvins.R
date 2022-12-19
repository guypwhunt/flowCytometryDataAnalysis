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

summary(combinedDf$X)

combinedDf <- combinedDf[grep("Rv", combinedDf$X),]
combinedDf$lipid <- combinedDf$X

combinedDf <- combinedDf[combinedDf$experiment %in%
                           combinedDf[combinedDf$adj.P.Val < 0.2 , "experiment"], ]


combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupFastBulbarEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowLimbLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowBulbarLate.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupFastLimbEarly.csv",]
combinedDf <- combinedDf[combinedDf$experiment != "FastVsSlowFastEarlyVsFastLate.csv",]
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
combinedDf <- combinedDf[!grepl("ControlVs", combinedDf$experiment, fixed = TRUE),]



combinedDf[combinedDf$experiment == "FastVsSlowFastLateVsFastEarly.csv", "experiment"] <- "Late Fast vs Early Fast"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowBulbarLate.csv", "experiment"] <- "Late Fast Bulbar vs Late Slow Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastLimbLate.csv", "experiment"] <- "Late Slow Limb vs Late Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowLimbLate.csv","experiment"] <- "Late Fast Limb vs Late Slow Limb"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbLateVsControl.csv","experiment"] <- "Late Limb vs Control"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbEarlyVsControl.csv","experiment"] <- "Early Limb vs Control"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbEarlyVsControl.csv","experiment"] <- "Early Limb vs Control"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlBulbarLateVsControl.csv","experiment"] <- "Late Bulbar vs Control"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlBulbarLateVsControl.csv","experiment"] <- "Early Bulbar vs Control"
combinedDf[combinedDf$experiment == "ProgressionVsControlSlowLateVsControl.csv","experiment"] <- "Late Slow vs Control"
combinedDf[combinedDf$experiment == "ProgressionVsControlSlowEarlyVsControl.csv","experiment"] <- "Early Slow vs Control"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastLateVsControl.csv","experiment"] <- "Late Fast vs Control"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastEarlyVsControl.csv","experiment"] <- "Early Fast vs Control"
combinedDf[combinedDf$experiment == "ProgressionVsControlFastEarlyVsControl.csv","experiment"] <- "Early Fast vs Control"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlFastLimbLateVsControl.csv","experiment"] <- "Late Fast Limb vs Control"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupSlowBulbarEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Slow Bulbar"
#combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastBulbarEarly.csv","experiment"] <- "Early Fast Limb vs Early Fast Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastLimbLate.csv","experiment"] <- "Early Fast Limb vs Late Fast Limb"



ggplot(combinedDf, aes(x = as.factor(lipid), y = -log10(adj.P.Val),
                      color = logFC)) +
  geom_point(alpha = 0.75, size=5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15)) +
  xlab("Specialized Pro-Resolving Mediator") + ylab("-log10(Adjusted P-Value)") +
  ylim(0, 3.5) + guides(color=guide_legend(title="Comparison")) +
  geom_hline(yintercept=0 - log10(0.05), linetype="dashed") +
  facet_wrap(~experiment) +
  guides(color = guide_colourbar(title="log2(Fold Change)")) +
  scale_colour_viridis_c()
