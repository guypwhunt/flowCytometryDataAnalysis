###############################
### LOAD REQUIRED PACKAGES ####
###############################

library(flowCore)
library(FlowSOM)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(RColorBrewer)
library(destiny)
library(uwot)
library(slingshot)
library(cytofkit2)
library(ggrepel)
library(tidyverse)


############################
######## LOAD DATA #########
############################

### Load the (transformed, normalized, unclustered) FCS files from the 'CSV_to_transformed_normalized_FCS' script
### Or load fcs files which where clustered in Cytosplore (in this case each fcs file is 1 cluster)

try(setwd("./flowCytometryData"))

dir.create("clusteringOutput", showWarnings = FALSE)

## Provide the directory of the fcs files
dirFCS <- paste0(getwd(), "/dataPPOutput")

## Optional: when loading clustered fcs files from cytosplore, provide the directory of the text file 'CSPLR_ST.txt'. Cytosplore exports this file upon running the HSNE. This file contains the decoding of the sample numbers.
pathST <- "X:/Users/guypw/OneDrive/Documents/flowCytometryDataAnalysis/HSNE_clusters_CD4/CSPLR_ST.txt"

## Defining a function to read multiple fcs files from a directory 'dir' into a single data.frame:
# NB: The column in the output named 'fileName' tracks the original file where each cell came from.
# Optionally perform remapping of column 'CSPLR_ST' holding cytosplore sample numbers to actual names:
read.flowdat <- function(dir,path_CSPLR_ST=""){
  # Read:
  filepaths <- list.files(path=dir,pattern = ".fcs", full.names=TRUE)
  flowset <- read.flowSet(files=filepaths, transformation=FALSE, truncate_max_range = FALSE)
  # Transform to data frame:
  x <- as.data.frame(exprs(as(flowset,'flowFrame')),stringsAsFactors=FALSE)
  # Map column 'Original' to filename (in this case holding clusters of HSNE):
  filenames <- gsub("[.fcs]","",list.files(path=dir,pattern = ".fcs", full.names=FALSE))
  names(filenames) <- sort(unique(x$Original))
  x$fileName <- filenames[as.character(x$Original)]
  # Remove column 'Original':
  x <- x[,-which(colnames(x)=="Original")]
  # Optionally remap Cytosplore sample tags to original filename:
  if(file.exists(path_CSPLR_ST)){
    # Read:
    sampID <- gsub(".fcs","",basename(sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[1])))
    names(sampID) <- sapply(strsplit(readLines(path_CSPLR_ST),": "),function(x) x[2])
    x$sampleID <- sampID[as.character(x$CSPLR_ST)]
  }
  return(x)
}

## Read fcs files
# In our example we will read the data which were clustered in Cytosplore (each fcs file is 1 cluster)
df <- read.flowdat(dir=dirFCS[1],path_CSPLR_ST = pathST)
gc()
write.csv(df, 'clusteringOutput/rawDf.csv')
gc()

# df <- read.csv('clusteringOutput/rawDf.csv')

df <- df[,14:28]

clinicalData <- read.csv('metadata/metadata.csv')

df <- merge(df, clinicalData, by.x = "fileName", by.y = "ï..id")

df["caseControl"][df["caseControl"] == "Case"] <- 1
df["caseControl"][df["caseControl"] == "Control"] <- 0

df["fastSlow"][df["fastSlow"] == "Fast"] <- 2
df["fastSlow"][df["fastSlow"] == "Slow"] <- 1
df["fastSlow"][df["fastSlow"] == "N/A"] <- 0
gc()

# http://www.sthda.com/english/articles/40-regression-analysis/168-multiple-linear-regression-in-r/

caseModel <- lm(caseControl ~ FSC.A + SSC.A + GPR32...AF488.A + FoxP3.PE.A +
              CD19...PE.CF595.A + IgD...PerCP.Cy5.5.A + IFNg.PE.Cy7.A +
              FPRL1...AF647.A + Zombie.NIR.A + IL.17...BV421.A +
              CD24...BV605.A + CD27...BV650.A, data = df)

print(summary(caseModel))

caseModelCoefficient <- as.data.frame(summary(caseModel)$coefficient)
caseModelCoefficient[order(caseModelCoefficient["Pr(>|t|)"]),]

print(confint(caseModel))

fastModel <- lm(fastSlow ~ FSC.A + SSC.A + GPR32...AF488.A + FoxP3.PE.A +
                  CD19...PE.CF595.A + IgD...PerCP.Cy5.5.A + IFNg.PE.Cy7.A +
                  FPRL1...AF647.A + Zombie.NIR.A + IL.17...BV421.A +
                  CD24...BV605.A + CD27...BV650.A, data = df)

print(summary(fastModel))

fastModelCoefficient <- as.data.frame(summary(fastModel)$coefficient)
fastModelCoefficient[order(fastModelCoefficient["Pr(>|t|)"]),]

print(confint(fastModel))
