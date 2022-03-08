library(R.utils)

workingDirectory <- getwd()

dataDirectorys <- c("/data/bCells","/data/monocytes",
                    "/data/senescence","/data/tCells")

for (directory in dataDirectorys) {
  try(setwd(paste0(workingDirectory,directory)))

  print(getwd())

  filenames <- list.files(pattern = ".csv")

  for (filename in filenames) {
    gzip(filename)
  }
}

setwd(workingDirectory)

