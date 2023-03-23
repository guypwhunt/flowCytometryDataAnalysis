try(source("R/01_functions.R"))

loadlibraries()

gpr18df <- read.csv("./data/isotypesGPR18/pairwiseTTestResults/senescence_ISO.csv")
gpr32df <- read.csv("./data/isotypesGPR32/rawSummaryStats/rawSummaryStats.csv")

df <- rbind(gpr18df, gpr32df)

df <- df[df$Marker %in% c("GPR18","Chem23", "GPR32...AF488.A", "FPRL1...AF647.A", "GPR32.AF488.A", "FPRL1.AF647.A"),]

df[df$Marker == "GPR32...AF488.A", "Marker"] <- "GPR32"
df[df$Marker == "GPR32.AF488.A", "Marker"] <- "GPR32"
df[df$Marker == "FPRL1...AF647.A", "Marker"] <- "FPRL1"
df[df$Marker == "FPRL1.AF647.A", "Marker"] <- "FPRL1"

df <- df[grepl("ISO", df$Type),]

df$colour <- NA

for(type in unique(df$Type)) {
  if (type == "bCells_ISO") {
    df[df$Type == type, "Type"] <- "B Cells"
  } else if (type == "monocytes_ISO") {
    df[df$Type == type, "Type"] <- "Monocytes"
  } else if (type == "tCells_ISO") {
    df[df$Type == type, "Type"] <- "T Cells"
  } else if (type == "senescence_ISO") {
    df[df$Type == type, "Type"] <- "Senescent T Cells"
  }
}


for(marker in unique(df$Marker)) {
  for (type in unique(df$Type)) {
    subselectedDf <-
      df[df$Marker == marker & df$Type == type, "Fold_Change"]
    try({
      if (min(quantile(subselectedDf, prob = c(.25))) > 0) {
        df[df$Marker == marker & df$Type == type, "colour"] <- "Good"
      } else {
        df[df$Marker == marker & df$Type == type, "colour"] <- "A bit dodgy"
      }
    })
  }
}

ggplot(df, aes(x=Type, y=Fold_Change, fill=colour)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Marker) +
  labs(y = "Fold Change", x = "Panel") +
  guides(fill=guide_legend(title="Signal")) +
  geom_hline(yintercept=0, linetype="dashed") #+
  # ylim(-2,2)
