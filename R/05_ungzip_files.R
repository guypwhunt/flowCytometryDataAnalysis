library(R.utils)

try(setwd("./flowCytometryData"))

gzFilenames <- list.files(pattern = ".gz")

for (gzFilename in gzFilenames) {
  gunzip(gzFilename)
}
