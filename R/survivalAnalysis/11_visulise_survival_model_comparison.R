library(dplyr)
library(ggplot2)

df <- fread(paste0("data/pValueAdjustmentsResults/", "diseaseDuration.csv"))

df <- as.data.frame(df)

df[1, ncol(df)] <- 1

# df <- na.omit(df)

colnames(df) <- c("Model", "LogLikelihood", "ChiSquared", "DegreesOfFreedom", "GlobalSignificance")

df[, "minusLogGlobalSignificance"] <- 0-log10(df[, "GlobalSignificance"])

df$Model <- factor(df$Model, levels = unique(df$Model))

ggplot(data=df, aes(x=minusLogGlobalSignificance, y=Model, fill = ChiSquared)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("-log10(P-Value)") +
  ylab("Survival Model") +
  guides(fill = guide_colourbar(title = "Chi-Squared"))


