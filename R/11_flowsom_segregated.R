set.seed(42)
library(flowCore)
library(FlowSOM)

# Read data
directory <- "data/bCells/dataPPOutput"

fileNames <- list.files(path = "data/bCells/dataPPOutput", pattern = ".fcs")

fileDirectories <- paste0(directory, "/", fileNames)

fSOM <- ReadInput(directory, compensate = FALSE, transform = FALSE, toTransform = c(), scale = FALSE)

#fSOM <- ReadInput(fileDirectories, compensate = FALSE, transform = FALSE, toTransform = c(), scale = FALSE)

# Building the self-organizing map
colnames(fSOM$data)

fSOM <- BuildSOM(fSOM, colsToUse = c(18, 23:24))

# Building the minimal spanning tree
fSOM <- BuildMST(fSOM)
PlotStars(fSOM)
PlotStars(fSOM, view = "grid")

tsne <- Rtsne::Rtsne(fSOM$map$codes, perplexity = 6)
PlotStars(fSOM, view = tsne$Y, maxNodeSize = 2)

# Plot marker
p <- PlotMarker(fSOM, "IgD...PerCP.Cy5.5.A")
print(p)

p <- PlotMarker(fSOM, "CD24...BV605.A")
print(p)

p <- PlotMarker(fSOM, "CD27...BV650.A")
print(p)

p <- PlotMarker(fSOM, "GPR32...AF488.A")
print(p)

p <- PlotMarker(fSOM, "FPRL1...AF647.A")
print(p)

# Plot Node number
PlotNumbers(fSOM)

# Compare markers across clusters
plot <- Plot2DScatters(fSOM,
                       channelpairs = list(c("CD27...BV650.A", "CD24...BV605.A")),
                       clusters = list(c(81, 82, 91, 92, 93)),
                       plotFile = NULL)
print(plot[[1]])

### Meta-clustering the data ###
metaClustering <- as.character(metaClustering_consensus(fSOM$map$codes,k = 7))

PlotLabels(fSOM, labels = metaClustering)

metaClustering_perCell <- GetMetaclusters(fSOM, metaClustering)
table(metaClustering_perCell)
metaClustering_perCell


# Detecting nodes with specific markers
query <- c("CD27...BV650.A" = "low",
           "CD24...BV605.A" = "high",
           "IgD...PerCP.Cy5.5.A" = "low")

query_res <- QueryStarPlot(fSOM, query, equalNodeSize = TRUE, plot = FALSE)
cellTypes <- factor(rep("Unlabeled", fSOM$map$nNodes),
                      levels=c("Unlabeled", "Immature B Cells"))
cellTypes[query_res$selected] <- "Immature B Cells"

p <- PlotStars(fSOM, backgroundValues=cellTypes,
               backgroundColor=c("#FFFFFF00","#0000FF"))
print(p)


#######
x <- readRDS("data/bCells/clusteringOutput/flowSom.rds")

query <- c("CD27...BV650.A" = "low",
           "CD24...BV605.A" = "high",
           "IgD...PerCP.Cy5.5.A" = "low")

query_res <- QueryStarPlot(x, query, equalNodeSize = TRUE, plot = FALSE)
cellTypes <- factor(rep("Unlabeled", x$map$nNodes),
                    levels=c("Unlabeled", "Immature B Cells"))
cellTypes[query_res$selected] <- "Immature B Cells"

p <- PlotStars(x, backgroundValues=cellTypes,
               backgroundColor=c("#FFFFFF00","#0000FF"))
print(p)
