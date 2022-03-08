###############################
### LOAD REQUIRED PACKAGES ####
###############################

library(flowCore)
library(Biobase)
library(dplyr)
library(flowVS)
library(flowStats)
library(R.utils)

###############################
##### LOAD PREPARED DATA ######
###############################

## Save the the compensated events in the gate of interest per individual sample as a csv file (for FlowJo users: select scale values).
# Set working directory:
setwd("./Senescence CSV")
rm(list = ls())

# Find file names of .csv files in the current working directory:
filenames <- list.files(pattern = ".csv")

# Verify:
filenames

columnNames <- c()

for (filename in filenames) {
  df <- read.csv(filename, sep=",", header=TRUE)
  print(filename)
  print(identical(colnames(df), columnNames))
  columnNames <- colnames(df)
}

rm(df)


## Defining a function to read a flow cytrometry file in csv format:
# Each row is a cell, each column is a parameter. In our experience, the flow cytometers sometimes output duplicate entries (listing the same cell twice), we remove these and report.
# Please check how your csv file is separated and adjust the sep argument in the function if necessary. In this example we import a semicolon separated file.
read.flow_csv <- function(pathIN){
  raw <- read.csv(pathIN, sep=",", header=TRUE, stringsAsFactors=FALSE)
  IND <- which(duplicated(raw))
  # Check for duplicates and report if found:
  if(any(duplicated(raw))){
    cat(paste0("=== Duplicate entries removed in [",pathIN,"]: ",length(IND)," ===\n"))
    print(head(raw[IND,]))
    cat("----\n")
  }
  return(unique(raw))
}

# Read all:
dfs <- sapply(filenames,read.flow_csv,simplify=FALSE)
gc()

##############################
#REWRITE TO FLOWFRAME/FLOWSET#
##############################

## Defining a function to rewrite a csv into a flowframe:
csv_2_ff <- function(dat){
  # Compute required metadata - column names with description - ranges, min, and max settings
  meta <- data.frame(name=dimnames(dat)[[2]],
                     desc=paste(dimnames(dat)[[2]]),
                     range =(apply(apply(dat,2,range),2,diff)),
                     minRange = apply(dat,2,min),
                     maxRange = apply(dat,2,max))
  # Create flowframe
  flowframef <- new("flowFrame",exprs=as.matrix(dat),parameters=AnnotatedDataFrame(meta))
  return(flowframef)
}

# rewrite to flowframe
dfs_ff = sapply(dfs,function(x) csv_2_ff(x),simplify=FALSE)
gc()
rm(dfs)

# rewrite to flowset
dfs_fs <- as(dfs_ff,"flowSet")
gc()
rm(dfs_ff)


###############################
####### TRANSFORMATION ########
###############################

## Each parameter of interest needs to be arcsinh transformed with an individual cofactor. The cofactor can be deduced from the size of the linear region around zero on a biexponential scale, as plotted in a histogram (in conventional gating software).
# Choose manual transformation or automated transformation (we prefer manual)
# Define parameters and cofactors for transformations:
#manualcofactors <- c(FSC.A=14000,SSC.A=24500,GPR32...AF488.A=170,FoxP3.PE.A=430,CD19...PE.CF595.A=1350,IgD...PerCP.Cy5.5.A=2.75,IFNg.PE.Cy7.A=24,FPRL1...AF647.A=2650,Zombie.NIR.A=1390,IL.17...BV421.A=115,
#                     CD24...BV605.A=1900,CD27...BV650.A=6.5)
#dfs_fs_t_manual <- transFlowVS(dfs_fs, channels=names(manualcofactors), cofactor=manualcofactors)
gc()

columnNames <- columnNames[c(16:(length(columnNames)-2))]
gc()
jpeg(file = "test.jpeg")
automatedcofactors <- estParamFlowVS(dfs_fs, columnNames) #this may take a while.
dev.off()
gc()

#auto
dfs_fs_t_auto <- transFlowVS(dfs_fs, channels=columnNames, cofactor=automatedcofactors)
gc()
rm(dfs_fs)
rm(automatedcofactors)


##############################
######## NORMALIZATION #######
##############################

## To correct for technical inter-sample variation we apply normalization by fdaNorm (which automatically detects the number of peaks)
# We continue with the manual transformed dataset
# Select the markers which require normalization (based on the densityplots you generated above). Be aware that you don't remove biological variation!
gc()
dfs_fs_t_auto_normfda <- warpSet(dfs_fs_t_auto, stains=columnNames)
gc()

##############################
####### EXPORT TO FCS ########
##############################

## The flowset (dfs_fs_t_auto_normfda) can be exported to individual fcs files
# Create an 'output' folder
gc()
dir.create("dataPPOutput", showWarnings = FALSE)
gc()

#Save flowframes wihtin flowset as fcs files using the flowCore package
write.flowSet(dfs_fs_t_auto_normfda, outdir='dataPPOutput', filename = paste0(gsub(".csv", ".fcs", sampleNames(dfs_fs_t_auto_normfda))))
gc()

##############################
########### Plots ############
##############################

# Pre-Normalized Plots
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

dir.create("figures", showWarnings = FALSE)
figureDirectory <- paste0(getwd(),"/figures/")

gc()
jpeg(file = paste0(figureDirectory,"densityPlotGPR32.jpeg"))
densityplot(~GPR32.AF488.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotKLRG1.jpeg"))
densityplot(~KLRG1.PE.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD4.jpeg"))
densityplot(~CD4.PE.CF594.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD8PerCP.jpeg"))
densityplot(~CD8.PerCP.Cy5.5.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCCR7.jpeg"))
densityplot(~CCR7.PE.Cy7.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotFPRL1.jpeg"))
densityplot(~FPRL1.AF647.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotZombie.jpeg"))
densityplot(~Zombie.NIR.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD27.jpeg"))
densityplot(~CD27.BV421.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotx.jpeg"))
densityplot(~x.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotPacific.jpeg"))
densityplot(~Pacific.Orange.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD45RA.jpeg"))
densityplot(~CD45RA.BV605.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD8BV650.jpeg"))
densityplot(~CD8.BV650.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"densityPlotCD28.jpeg"))
densityplot(~CD28.BV785.A, dfs_fs_t_auto, main="auto")
dev.off()
gc()

# Normalized Plots
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotGPR32.jpeg"))
densityplot(~GPR32.AF488.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotKLRG1.jpeg"))
densityplot(~KLRG1.PE.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD4.jpeg"))
densityplot(~CD4.PE.CF594.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD8PerCP.jpeg"))
densityplot(~CD8.PerCP.Cy5.5.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCCR7.jpeg"))
densityplot(~CCR7.PE.Cy7.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotFPRL1.jpeg"))
densityplot(~FPRL1.AF647.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotZombie.jpeg"))
densityplot(~Zombie.NIR.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD27.jpeg"))
densityplot(~CD27.BV421.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotx.jpeg"))
densityplot(~x.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotPacific.jpeg"))
densityplot(~Pacific.Orange.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD45RA.jpeg"))
densityplot(~CD45RA.BV605.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD8BV650.jpeg"))
densityplot(~CD8.BV650.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
jpeg(file = paste0(figureDirectory,"normalisedDensityPlotCD28.jpeg"))
densityplot(~CD28.BV785.A, dfs_fs_t_auto_normfda, main="auto")
dev.off()
gc()
