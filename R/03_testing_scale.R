
data = c(1200,34567,3456,12,3456,0985,1211)
summary(data)
library(caret)
process <- preProcess(as.data.frame(data), method=c("range"))

norm_scale <- predict(process, as.data.frame(data))

#######
df <- read.csv("C:/Users/guypw/OneDrive/Documents/resolvinAnalysis/data/bCells/clusteringOutput/flowSomDf.csv")

colnames(df)

columnNames <- c("GPR32...AF488.A", "FPRL1...AF647.A", "IgD...PerCP.Cy5.5.A", "CD24...BV605.A","CD27...BV650.A")

df <- df[, columnNames]

minMaxDF <- df
scaleDF <- df
scaledDF <- scale(df)


for (col in columnNames) {
  dfMax <- max(minMaxDF[, col])
  dfMin <- min(minMaxDF[, col])

  scaledColumn <- (minMaxDF[, col] - dfMin) /(dfMax-dfMin)

  d <- density(scaledColumn)
  print(plot(d))

  minMaxDF[, col] <- scaledColumn
}
