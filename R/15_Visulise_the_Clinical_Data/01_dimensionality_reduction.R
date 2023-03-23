# install.packages("umap")
library(umap)
library(readxl)

sampleInformation <- read_excel("C:/Users/guypw/OneDrive/Documents/resolvinAnalysis/data/metadata/clinicalData.xlsx")
sampleInformation <- as.data.frame(sampleInformation)
sampleInformation <- sampleInformation[, -1]

sampleInformation <- sampleInformation[sampleInformation$visit == 1,]
# sampleInformation <- sampleInformation[sampleInformation$outcome == "dead",]
colnames(sampleInformation)

sampleInformation <- sampleInformation[,c("fastSlow", "fastSlowFrs", "diagnosticDelayInYears", "diseaseDurationInYears", "delataAlsfrsRScore", "ageAtOnset", "BulbarLimb"
                                          )]
sampleInformation <- na.omit(sampleInformation)
sampleInformation$BulbarLimb <- as.numeric(as.factor(sampleInformation$BulbarLimb))


umapResults <- umap(sampleInformation[, c("diagnosticDelayInYears", "diseaseDurationInYears", "delataAlsfrsRScore")])

umapResults$layout <- as.data.frame(umapResults$layout)
colnames(umapResults$layout) <- c("UMAP1", "UMAP2")
umapResults$layout$Progression <- sampleInformation[, "fastSlow"]

# x <- sampleInformation[umapResults$layout$UMAP1 >0,]
# summary(x[x$fastSlow == "Fast","diseaseDurationInYears"])

ggplot(umapResults$layout, aes(UMAP1, UMAP2, colour = Progression)) +
  geom_point()

umapResults$layout$Progression <- sampleInformation[, "fastSlowFrs"]

ggplot(umapResults$layout, aes(UMAP1, UMAP2, colour = Progression)) +
  geom_point()

