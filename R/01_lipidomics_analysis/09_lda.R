library(limma)
library(umap)
library(readxl)
library(data.table)
library(dplyr)
library(ggrepel)
library(tibble)
library(factoextra)
library(tidyr)
library(ltm)
library(MASS)
library(tidyverse)
library(caret)

try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

# metacluster_colours <-
#   as.vector(colortools::wheel('#8b0000', num = 4))

visit <- 1

sampleInformation <-
  read_excel("data/metadata/clinicalData.xlsx") %>%
  as.data.frame() %>%
  dplyr::filter(visit == !!visit | caseControl == "Control") %>%
  dplyr::filter(experiment == "lipidomics") %>%
  dplyr::select(classification, fastSlow
  )

columnName <- "fastSlow"

factorLevels <- c("Baseline A-F", "Baseline A-S", "NCC")

resolvinsOnly <- FALSE

performLDA(sampleInformation, columnName, visit, factorLevels, resolvinsOnly)

############
visit <- 2

sampleInformation <-
  read_excel("data/metadata/clinicalData.xlsx") %>%
  as.data.frame() %>%
  dplyr::filter(visit == !!visit | caseControl == "Control") %>%
  dplyr::filter(experiment == "lipidomics") %>%
  dplyr::select(classification, fastSlow
  )

columnName <- "fastSlow"

factorLevels <- c("V2 A-F", "V2 A-S", "NCC")

performLDA(sampleInformation, columnName, visit, factorLevels, resolvinsOnly)


###############
visit <- 12

sampleInformation <-
  read_excel("data/metadata/clinicalData.xlsx") %>%
  as.data.frame() %>%
  dplyr::filter(experiment == "lipidomics") %>%
  dplyr::select(classification, earlyLate
  )

columnName <- "earlyLate"

factorLevels <- c("V2 pwALS", "Baseline pwALS", "NCC")

performLDA(sampleInformation, columnName, visit, factorLevels, resolvinsOnly)
