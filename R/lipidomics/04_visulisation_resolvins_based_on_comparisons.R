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

summary(combinedDf$logFC)

df <- combinedDf

combinedDf$lipid <- factor(combinedDf$X, levels = unique(combinedDf[order(combinedDf$X), "X"]))

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


combinedDf <- combinedDf[!grepl("ControlVs", combinedDf$experiment, fixed = TRUE),]



combinedDf[combinedDf$experiment == "FastVsSlowFastEarlyVsFastLate.csv", "experiment"] <- "Early Fast vs Late Fast"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupSlowBulbarLate.csv", "experiment"] <- "Late Fast Bulbar vs Late Slow Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowLimbLateVsgroupFastLimbLate.csv", "experiment"] <- "Late Slow Limb vs Late Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbLateVsgroupSlowLimbLate.csv","experiment"] <- "Late Fast Limb vs Late Slow Limb"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlLimbLateVsControl.csv","experiment"] <- "Late Limb vs Control"
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
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastBulbarEarly.csv","experiment"] <- "Early Fast Limb vs Early Fast Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastLimbEarly.csv","experiment"] <- "Early Fast Bulbar vs Early Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupFastLimbLate.csv","experiment"] <- "Early Fast Limb vs Late Fast Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowBulbarEarlyVsControl.csv","experiment"] <- "Early Slow Bulbar vs Control"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarEarlyVsBulbarLate.csv","experiment"] <- "Early Bulbar vs Late Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupSlowBulbarLate.csv","experiment"] <- "Early Slow Bulbar vs Late Slow Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowLimbLate.csv" ,"experiment"] <- "Late Slow Bulbar vs Late Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarLateVsgroupSlowLimbLate.csv" ,"experiment"] <- "Late Slow Bulbar vs Late Slow Limb"
combinedDf[combinedDf$experiment == "FastVsSlowFastEarlyVsSlowEarly.csv" ,"experiment"] <- "Early Fast vs Early Slow"
combinedDf[combinedDf$experiment == "FastVsSlowFastLateVsSlowLate.csv" ,"experiment"] <- "Late Fast vs Late Slow"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowBulbarEarlyVsgroupSlowLimbEarly.csv" ,"experiment"] <- "Early Slow Bulbar vs Early Slow Limb"
combinedDf[combinedDf$experiment == "SightOfOnsetVsControlBulbarEarlyVsControl.csv" ,"experiment"] <- "Early Bulbar vs Control"
combinedDf[combinedDf$experiment == "CaseVsControlEarlyVsControl.csv" ,"experiment"] <- "Early vs Control"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarEarlyVsgroupFastBulbarLate.csv" ,"experiment"] <- "Early Fast Bulbar vs Late Fast Bulbar"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastBulbarLateVsgroupFastLimbLate.csv" ,"experiment"] <- "Late Fast Bulbar vs Late Fast Limb"
combinedDf[combinedDf$experiment == "EarlyVsLateEarlyVsLate.csv" ,"experiment"] <- "Early vs Late"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupSlowLimbEarlyVsgroupSlowLimbLate.csv" ,"experiment"] <- "Early Slow Limb vs Late Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetgroupFastLimbEarlyVsgroupSlowLimbEarly.csv" ,"experiment"] <- "Early Fast Limb vs Early Slow Limb"
combinedDf[combinedDf$experiment == "ProgressionAndSightOfOnsetVsControlSlowBulbarLateVsControl.csv" ,"experiment"] <- "Slow Bulbar vs Control"
combinedDf[combinedDf$experiment == "CaseVsControlLateVsControl.csv" ,"experiment"] <- "Late vs Control"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarEarlyVsLimbEarly.csv" ,"experiment"] <- "Early Bulbar vs Early Limb"
combinedDf[combinedDf$experiment == "BulbarVsLimbBulbarLateVsLimbLate.csv" ,"experiment"] <- "Late Bulbar vs Late Limb"
combinedDf[combinedDf$experiment == "BulbarVsLimbLimbEarlyVsLimbLate.csv" ,"experiment"] <- "Early Limb vs Late Limb"
combinedDf[combinedDf$experiment == "FastVsSlowSlowEarlyVsSlowLate.csv" ,"experiment"] <- "Early Slow vs Late Slow"

combinedDf$regulation <- "Not Statistically Different"

try({combinedDf[combinedDf$adj.P.Val < 0.05 &
                  combinedDf$logFC < 0, "regulation"] <-
  "Downregulated"})
try({combinedDf[combinedDf$adj.P.Val < 0.05 &
                  combinedDf$logFC > 0, "regulation"] <-
  "Upregulated"})

scaleLimits <-
  ceiling(max(abs(min(
    na.omit(df$logFC)
  )), max(
    na.omit(df$logFC)
  )))

regulation <- c(
  "Upregulated"  = "Red",
  "Not Statistically Different" = "Grey",
  "Downregulated" = "Blue"
)

combinedDf$regulation <- factor(combinedDf$regulation, levels = names(regulation))

# Early vs Controls
comparisons1 <- c(
  "Early vs Control" = 21,
  "Early Fast vs Control"  = 22,
  "Early Slow vs Control" = 23,
  "Early Bulbar vs Control"  = 24,
  "Early Limb vs Control" = 25
)

# Late vs Controls
comparisons2 <- c(
  "Late vs Control" = 21,
  "Late Fast vs Control"  = 22,
  "Late Slow vs Control" = 23,
  "Late Bulbar vs Control"  = 24,
  "Late Limb vs Control" = 25
)

# ALS Progression Subgroup Comparisons
comparisons3 <- c(
  "Early Fast vs Early Slow" = 21,
  "Early Fast vs Late Fast"  = 22,
  "Early Slow vs Late Slow" = 23,
  "Late Fast vs Late Slow" = 24
)

# ALS Site of Onset Subgroup Comparisons
comparisons4 <- c(
  "Early Bulbar vs Early Limb" = 21,
  "Early Bulbar vs Late Bulbar"  = 22,
  "Early Limb vs Late Limb" = 23,
  "Late Bulbar vs Late Limb" = 24
)

comparisonsList <- list(comparisons1, comparisons2, comparisons3, comparisons4)

combinedDfGoldenSource <- combinedDf

for(comparisons in comparisonsList) {
  combinedDf <- combinedDfGoldenSource
  combinedDf <-
    combinedDf[combinedDf$experiment %in% names(comparisons),]
  combinedDf <- combinedDf[grep("Rv", combinedDf$X), ]
  combinedDf <- combinedDf[nchar(combinedDf$X) == 4, ]
  combinedDf$experiment <-
    factor(combinedDf$experiment, levels = names(comparisons))

  print(
    ggplot(
      combinedDf,
      aes(
        x = lipid,
        y = -log10(P.Value),
        fill = logFC,
        color = regulation,
        shape = experiment
      )
    ) +
      geom_point(
        alpha = 1,
        size = 5,
        stroke = 1,
        position = position_jitter(width = 0.1, height = 0)
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1,
          size = 8
        ),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.justification = "top"
      ) +
      scale_shape_manual(values = comparisons) +
      scale_colour_manual(values = regulation, drop = FALSE) +
      xlab("Resolvins") + ylab("-log10(P-Value)") +
      guides(color = guide_legend(title = "Comparison")) +
      geom_hline(yintercept = 0 - log10(0.05), linetype = "dashed") +
      guides(
        fill = guide_colourbar(title = "log2(Fold Change)", order = 2),
        color = guide_legend(
          title = "Expression",
          order = 1,
          override.aes = list(shape = 1, size = 4)
        ),
        shape = guide_legend(
          title = "Comparison",
          order = 3,
          override.aes = list(size = 4)
        )
      ) +
      ylim(0, ceiling(max(
        0 - log10(df$P.Value)
      ))) +
      scale_fill_gradientn(
        limits = c(0 - scaleLimits, scaleLimits),
        colours = c(
          "#0000FF",
          "#2E2EFF",
          "#5C5CFF",
          "#8A8AFF",
          "#ffffff",
          "#FF8A8A",
          "#FF5C5C",
          "#FF2E2E",
          "#FF0000"
        )
      )
  )
}
