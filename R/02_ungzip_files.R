library(R.utils)

workingDirectory <- getwd()

dataDirectorys <- c("/data/bCells","/data/monocytes",
                    "/data/senescence","/data/tCells")

for (directory in dataDirectorys) {
  try(setwd(paste0(workingDirectory,directory)))

  gzFilenames <- list.files(pattern = ".gz")

  for (gzFilename in gzFilenames) {
    gunzip(gzFilename)
  }
}

setwd(workingDirectory)
