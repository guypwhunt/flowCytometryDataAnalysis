## INSTALL REQUIRED PACKAGES ##
###############################
install.packages(c("tidyverse","stringr","stringi","knitr","roxygen2","BiocManager","dplyr","R.utils","reshape2","ggplot2","uwot","ggrepel","dplyr","ggplot2","scales","reshape2","RColorBrewer","devtools"))
install.packages(c("locfit", "hdrcde", "rainbow", "fds", "fda", "flowStats", "openCyto", "CytoML"))

BiocManager::install("flowStats")
BiocManager::install("openCyto")
BiocManager::install("CytoML")
BiocManager::install("Biobase")
BiocManager::install("flowCore")
BiocManager::install("flowVS")
BiocManager::install("flowStats")
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
