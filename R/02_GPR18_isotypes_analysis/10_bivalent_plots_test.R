try(source("R/01_functions.R"))

loadlibraries()

workingDirectory <- getwd()

directoryName <- "bCells_ISO"

setwd(paste0("./data/isotypes/", directoryName))

# Create an 'output' folder
gc()
dir.create("results", showWarnings = FALSE)
dir.create("transformedData", showWarnings = FALSE)

filenames <- list.files(pattern = ".csv")

df <- read.csv(filenames[1])

q <- quantile(df[, "GPR32...AF488.A"], c(0.01, 0.99))

df <- df[df[, "GPR32...AF488.A"] < max(q) &
           df[, "GPR32...AF488.A"] > min(q), ]

d <- density(df[, "GPR32...AF488.A"])

min(d$x)
plot(d)

colnames(df)
df[df[, "GPR32...AF488.A"] >=10, "GPR32...AF488.A"]  <- log(df[df[, "GPR32...AF488.A"] >=10, "GPR32...AF488.A"]) + 10

d2 <- density(df[, "GPR32...AF488.A"])
min(d2$x)

plot(d2)
