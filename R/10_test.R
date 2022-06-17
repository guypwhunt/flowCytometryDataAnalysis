flowSom <- readRDS("data/bCells/clusteringOutput/flowSom.rds")

attributes(flowSom)
df <- read.csv("data/bCells/clusteringOutput/fastPGDf.csv", nrow = 10)
colnames(df)
#flowSom$pattern
#flowSom$compensate
#flowSom$spillover
#flowSom$transform
#flowSom$toTransform
#flowSom$transformFunction
#flowSom$transformList
#flowSom$scale
#flowSom$info


### Maybe I could update the column names before visualization
flowSom$prettyColnames
#flowSom$data
### Unknown
flowSom$metaData
### Maybe some interesting stats
flowSom$map
attributes(flowSom$map)
# outliers
flowSom$outliers
# MST
flowSom$MST
# metaclustering
flowSom$metaclustering

PlotPies(flowSom, flowSom$map$mapping)


df <- read.csv("data/bCells/dataPPOutput/columnsOfInterestDf.csv")
colnames(df)

df2 <- df[, c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A","CD27...BV650.A")]

set.seed(31)

tmp <- NULL
for (k in 1:10){
  tmp[k] <- kmeans(df2, k, nstart = 30)
}
df3 <- data.frame(tmp)
# add a prefix to the column names
colnames(df3) <- seq(1:10)
colnames(df3) <- paste0("k",colnames(df3))
# get individual PCA
df.pca <- prcomp(df3, center = TRUE, scale. = FALSE)
ind.coord <- df.pca$x
ind.coord <- ind.coord[,1:2]
df3 <- bind_cols(as.data.frame(df3), as.data.frame(ind.coord))

library(clustree)
clustree(df3, prefix = "k")

# https://www.r-bloggers.com/2019/01/10-tips-for-choosing-the-optimal-number-of-clusters/#:~:text=The%20NbClust%20package%20provides%2030,distance%20measures%2C%20and%20clustering%20methods.&text=This%20suggest%20the%20optimal%20number%20of%20clusters%20is%203.
