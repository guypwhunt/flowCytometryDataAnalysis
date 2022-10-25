try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

library(fpc)

directoryName <- "bCells"
columnNames <- bCellsClusteringColumnNames

df <- read.csv(paste0("./data/", directoryName, '/clusteringOutput/clusteringOutputs.csv'))

dirFCS <- paste0("./data/", directoryName, "/dataPPOutput/scaledFcs")

fcsDf <- read.twoflowdat(dir = dirFCS[1], path_CSPLR_ST = "")

cf1 <- clusterboot(df[, bCellsClusteringColumnNames],B=3,bootmethod=
                     c("boot","noise","jitter"),clustermethod=kmeansCBI,
                   krange=5,seed=15555)

print(cf1)
plot(cf1)



require(graphics)

face <- rFace(50,dMoNo=2,dNoEy=0,p=2)

cf2 <- clusterboot(dist(as.data.frame(face)),B=3,bootmethod=
                     "subset",clustermethod=disthclustCBI,
                   k=5, cut="number", method="average", showplots=TRUE, seed=15555)
print(cf2)
