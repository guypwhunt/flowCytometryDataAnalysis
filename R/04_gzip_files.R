library(R.utils)

try(setwd("./flowCytometryData"))

filenames <- list.files(pattern = ".csv")

for (gzFilename in filenames) {
  gzip(gzFilename)
}
