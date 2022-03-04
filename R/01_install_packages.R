## INSTALL REQUIRED PACKAGES ##
###############################
install.packages("BiocManager")
install.packages("dplyr")
BiocManager::install("Biobase")
BiocManager::install("flowCore")
BiocManager::install("flowVS")
BiocManager::install("flowStats")

install.packages("reshape2")
install.packages("ggplot2")
install.packages("uwot")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("scales")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("BiocManager")

BiocManager::install("FlowSOM")
BiocManager::install("slingshot")
BiocManager::install("flowCore")
BiocManager::install("SingleCellExperiment")

library(devtools)
devtools::install_github("JinmiaoChenLab/cytofkit2")
devtools::install_github('flying-sheep/knn.covertree')
devtools::install_github('theislab/destiny')
devtools::install_github('sararselitsky/FastPG')


gc()
