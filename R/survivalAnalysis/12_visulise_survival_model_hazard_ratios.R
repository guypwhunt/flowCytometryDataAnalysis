library(dplyr)
library(ggplot2)

figureDirectory <- "data/pValueAdjustmentsResults/"

fileNames <- list.files(figureDirectory, pattern ="GPR")

fileNames <- append(fileNames, "DiseaseDurationClinicalModel.csv")

fileNamesPath <- paste0(figureDirectory, fileNames)

dfs <- lapply(fileNamesPath, read.csv)

names(dfs) <- fileNames

i <- 1

for (df in dfs) {
  df[, "ModelName"] <- names(dfs)[i]

  dfs[[i]] <- df

  i <- i + 1
}

df <- do.call(rbind, dfs)

df <- as.data.frame(df)

unique(df$ModelName)

df[df$ModelName == "DiseaseDurationClinicalModel.csv" , "ModelName"] <- "Clinical"
df[df$ModelName == "GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR32 (FlowSOM)"
df[df$ModelName == "GPR32clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR32 (Phenograph)"
df[df$ModelName == "GPR18meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 (FlowSOM)"
df[df$ModelName == "GPR18clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 (Phenograph)"
df[df$ModelName == "GPR18GPR32meta_clusters_flowsomDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 & GPR32 (FlowSOM)"
df[df$ModelName == "GPR18GPR32clusters_phenographDiseaseDurationBiologicalModel.csv" , "ModelName"] <- "GPR18 & GPR32 (Phenograph)"

df$ModelName <- factor(df$ModelName, levels = c("Clinical", "GPR18 (FlowSOM)", "GPR18 (Phenograph)", "GPR32 (FlowSOM)", "GPR32 (Phenograph)", "GPR18 & GPR32 (FlowSOM)", "GPR18 & GPR32 (Phenograph)"))

colnames(df)[6] <- "p.value"

df$lower <- df$Hazard.Ratio - df$Confidence.Interval
df$upper <- df$Hazard.Ratio + df$Confidence.Interval


dfPositives <- df
# dfPositives$Hazard.Ratio <- log10(dfPositives$Hazard.Ratio)
dfPositives$lower <- dfPositives$Regression.Coefficient - dfPositives$Confidence.Interval
dfPositives$upper <- dfPositives$Regression.Coefficient + dfPositives$Confidence.Interval

dfPositives <- dfPositives[dfPositives$ModelName %in% c("Clinical", "GPR32 (Phenograph)", "GPR18 (Phenograph)", "GPR18 & GPR32 (Phenograph)"), ]

bCells <- c(
  "F-B (Median GPR18)",
  "F-B (Median GPR32)",
  "CD24+US-B (Median GPR32)",
  "CD24-US-B (Median GPR18)"
)

monocytes <- c(
  "HLA-DR-IM-M (Median GPR32)",
  "CD11b Low HLA-DR-aCD11b+IM-M (Median GPR32)",
  "HLA-DR-C-M (Median GPR32)",
  "CD11b Low HLA-DR-aCD11b+C-M (Median GPR32)"
)

cells <- list(bCells, monocytes)

for (cell in cells) {
  dfPositivesSub <- dfPositives[dfPositives$Variable %in% cell,]
  dfPositivesSub$Variable <-
    factor(dfPositivesSub$Variable, levels = cell)

  print(
    ggplot(
      data = dfPositivesSub,
      aes(
        y = Variable,
        x = Hazard.Ratio,
        xmin = lower,
        xmax = upper
      )
    ) +
      geom_point(size = 3) +
      geom_errorbarh(height = .1) +
      facet_wrap( ~ ModelName) +
      ylab("Cell Populations") +
      xlab("Regression Coefficient [95% Confidence Interval]") +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      geom_vline(xintercept = 0, linetype="dotted") +
      scale_x_continuous(trans = "log", breaks = scales::breaks_extended())
  )
}


# fulldata <- as_tibble(df)
#
# fulldata <- fulldata %>%
#   mutate(lab=case_when(p.value<0.0001 ~ "****",
#                        p.value<0.001 ~ "***",
#                        p.value<0.01 ~ "**",
#                        p.value<0.05 ~ "*",
#                        p.value<0.1 ~ "+",
#                        TRUE ~ ''))
#
#
#
#
# # fulldata <- fulldata %>%
# #   filter(!(model %in% c("Alzheimer's Disease","Frontotemporal dementia")))
#
# total_plot <- ggplot(fulldata,aes(x=Variable,y=exp(Hazard.Ratio)#,color=design, shape=covars
#                                  ))+
#   geom_point(position=position_dodge(width = 0.6))+
#   geom_errorbar(aes(ymin=exp(Hazard.Ratio - Confidence.Interval*0.8),ymax=exp(Hazard.Ratio + Confidence.Interval)),width=0,
#                 position=position_dodge(width = 0.6))+
#   facet_wrap(~ModelName,ncol=4#,scales = "free"
#              )+
#   geom_text(aes(y=exp(max(df$Confidence.Interval) + max(df$Hazard.Ratio))* 0.5,label=lab),show.legend = FALSE) + #Only show significance for the simple model
#   theme_bw() +
#   geom_hline(yintercept = 1,lty=2) +
#   ylab("Odds ratio [95% confidence interval]") +
#   xlab("Variable")+
#   #labs(caption = "P: **** < 0.0001 < *** < 0.001 < ** < 0.01 < * < 0.05 < + < 0.1")+
#   theme(legend.position = c(1, 0),
#         legend.justification = c(1.60,-0.10),
#         plot.caption = element_text(hjust = 1,size=10)#,
#         #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
#   )
#
# print(total_plot)
#
# # ##### SAVE THE PLOT
# # ggsave(file.path(path,"Flippedtotalprs_fullnew.pdf"),total_plot,device="pdf",units="mm",width=250,height=150)
# # #ggsave("~/Downloads/Flippedtotalprs.pdf",total_plot,device="pdf",units="mm",width=150,height=175)
# #
