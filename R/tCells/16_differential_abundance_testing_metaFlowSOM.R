# https://bioconductor.org/packages/devel/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html
try(source("R/01_functions.R"))

loadlibraries()

# variables
directoryName <- "tCells"

columnNames <- c("fileName", "CD127.BV510.A", "CD8.BV650.A",
                      "CD25.BV786.A", "FoxP3.PE.A", "CD45RO.PE.CF595.A",
                      "CD4.PerCP.Cy5.5.A", "GPR32.AF488.A", "FPRL1.AF647.A")



clusterName <- "meta_clusters_flowsom"

performAllDifferentialAbundanceTests(directoryName, columnNames, clusterName)
