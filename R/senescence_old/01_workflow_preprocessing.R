try(source("R/01_functions.R"))

loadlibraries()

senescenceDirectoryName <- "senescence"
senescenceColumnNames <- c("CD27.BV421.A",
                           "CD45RA.BV605.A",
                           "CD28.BV785.A",
                           "KLRG1.PE.A",
                           "CD4.PE.CF594.A",
                           "CD8.PerCP.Cy5.5.A",
                           "CD8.BV650.A",
                           "CCR7.PE.Cy7.A",
                           "GPR32.AF488.A",
                           "FPRL1.AF647.A",
                           "Zombie.NIR.A"
                           )
automatedcofactors <- c(5975.868, 100, 41.4155, 50,
                        5000, 10, 141.8866, 25.12086,
                        45.30041, 101.9894, 105.3142)

automatedcofactors_originals <- c(5975.868290, 4937.436002, 41.415505, 0.361972, 8742.139346, 3728.266090, 141.886562,
  25.120859, 45.300406, 101.989383, 105.314192)

automatedcofactors <- c(5975.868290, 4937.436002, 41.415505, 0.361972, 8742.139346, 1500, 141.886562,
                                  25.120859, 45.300406, 101.989383, 105.314192)


test <- FALSE

directoryName <- senescenceDirectoryName
columnNames <- senescenceColumnNames
gateColumns <- data.frame (CD8.PerCP.Cy5.5.A = c(2),
                           CD4.PE.CF594.A = c(1.5))
gate <- TRUE
gateTogether <- TRUE

preprocessing(directoryName, columnNames, test, gate, gateTogether,
              gateColumns, automatedcofactors = automatedcofactors)
