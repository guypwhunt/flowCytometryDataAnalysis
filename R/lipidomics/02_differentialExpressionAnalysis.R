library(limma)
library(umap)
library(readxl)


adjustment <- "fdr"

ex <- read.csv(
  "data/lipidomics/normalisedasinhTransformedExpressionDataRawOutliersAndDuplicatesRemoved.csv",
  row.names = 1
)

clinical <-
  read_excel("data/metadata/clinicalData.xlsx")
clinical <- as.data.frame(clinical)
clinical <- clinical[clinical$experiment == "lipidomics", ]

ex <- ex[, order(colnames(ex))]
ex <- ex[complete.cases(ex), ]

clinical <-
  clinical[which(clinical$classification %in% colnames(ex)), ]
rownames(clinical) <- clinical$classification
clinical <- clinical[colnames(ex), ]

all(colnames(ex) == row.names(clinical))

performDifferentialExpression <-
  function(ex,
           design,
           cont.matrix,
           experimentName,
           clinical = NULL,
           pairwiseTest = FALSE) {
    # put this in a function

    # calculate precision weights and show plot of mean-variance trend
    v <- vooma(ex, design, plot = TRUE)

    if (pairwiseTest) {
      corfit <- duplicateCorrelation(v, design, block = clinical$patient_id)

      fit <-
        lmFit(v,
              design,
              block = clinical$patient_id,
              correlation = corfit$consensus)
    } else {
      fit  <- lmFit(v)
    }

    fit2 <- contrasts.fit(fit, cont.matrix)

    # compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)

    # Visualize and quality control test results.
    # Build histogram of P-values for all genes. Normal test
    # assumption is that most genes are not differentially expressed.
    figureDirectory <- "data/lipidomics/Results/"

    par(mar = c(1, 1, 1, 1))

    jpeg(file = paste0(figureDirectory,
                       "histogram",
                       experimentName,
                       ".jpeg"))

    dev.off()
    gc()

    # summarize test results as "up", "down" or "not expressed"
    dT <-
      decideTests(fit2, adjust.method = adjustment, p.value = 0.05)


    jpeg(file = paste0(figureDirectory,
                       "qqPlot",
                       experimentName,
                       ".jpeg"))

    # create Q-Q plot for t-statistic
    t.good <- which(!is.na(fit2$F)) # filter out bad probes
    qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")
    abline(0, 1)

    dev.off()


    # volcano plot (log P-value vs log fold change)
    colnames(fit2) # list contrast names
    ct <- 1        # choose contrast of interest

    for (n in seq(ncol(fit2))) {
      jpeg(file = paste0(
        figureDirectory,
        "volcanoPlot",
        colnames(fit2)[n],
        experimentName,
        ".jpeg"
      ))
      volcanoplot(
        fit2,
        coef = n,
        main = colnames(fit2)[n],
        pch = 20,
        highlight = length(which(dT[, n] != 0)),
        names = rep('+', nrow(fit2))
      )
      dev.off()

      jpeg(
        file = paste0(
          figureDirectory,
          "meanDifferencePlot",
          colnames(fit2)[n],
          experimentName,
          ".jpeg"
        )
      )
      # MD plot (log fold change vs mean log expression)
      # highlight statistically significant (p-adj < 0.05) probes
      plotMD(
        fit2,
        column = n,
        status = dT[, n],
        legend = F,
        pch = 20,
        cex = 1
      )
      abline(h = 0)
      dev.off()

      tT <- topTable(
        fit2,
        adjust = adjustment,
        sort.by = "P",
        number = Inf,
        coef = n
      )

      write.csv(tT,
                paste0(figureDirectory, experimentName, colnames(fit2)[n], ".csv"))

      significantResults <- tT[tT$adj.P.Val < 0.05, ]

      if (nrow(significantResults) > 0) {
        print(colnames(fit2)[n])
        print(significantResults)
      }
    }
  }

nonControlClinical <- clinical[clinical$caseControl != "Control",]
nonControlClinical <-
  nonControlClinical[which(nonControlClinical$patient_id %in% nonControlClinical[duplicated(nonControlClinical$patient_id), "patient_id"]),]

nonControlEx <- ex[, colnames(ex) %in% rownames(nonControlClinical)]
nonControlEx <- ex[, colnames(ex) %in% rownames(nonControlClinical)]


### Case vs Control Experiment
clinical$group <- clinical$caseControlExperiment
nonControlClinical$group <- nonControlClinical$caseControlExperiment
experimentName <- "CaseVsControl"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + batchControl
  ,
  clinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  ControlVsEarly = groupControl - groupEarly,
  ControlVsLate = groupControl - groupLate,
  EarlyVsControl = groupEarly - groupControl,
  LateVsControl = groupLate - groupControl,
  levels = design
)

performDifferentialExpression(ex, design, cont.matrix, experimentName)

# No Controls
experimentName <- "EarlyVsLate"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + timeFromOnsetToVisitInYears
  + timeFromVisit1InYears
  + BulbarLimb
  + fastSlow
  + batchControl
  ,
  nonControlClinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  EarlyVsLate = groupEarly - groupLate,
  LateVsEarly = groupLate - groupEarly,
  levels = design
)

performDifferentialExpression(nonControlEx,
                              design,
                              cont.matrix,
                              experimentName,
                              nonControlClinical,
                              TRUE)


### Progression Experiment
clinical$group <- clinical$progressionExperiment
nonControlClinical$group <- nonControlClinical$progressionExperiment

experimentName <- "ProgressionVsControl"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + batchControl
  ,
  clinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  ControlVsSlowLate =   groupControl - groupSlowLate,
  ControlVsSlowEarly = groupControl - groupSlowEarly,
  ControlVsFastEarly = groupControl - groupFastEarly,
  ControlVsFastLate = groupControl - groupFastLate,
  SlowLateVsControl = groupSlowLate - groupControl,
  FastEarlyVsControl = groupFastEarly - groupControl,
  SlowEarlyVsControl = groupSlowEarly - groupControl,
  FastLateVsControl = groupFastLate - groupControl,
  levels = design
)

performDifferentialExpression(ex, design, cont.matrix, experimentName)

experimentName <- "FastVsSlow"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + timeFromOnsetToVisitInYears
  + timeFromVisit1InYears
  + BulbarLimb
  + batchControl
  ,
  nonControlClinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  SlowLateVsFastLate = groupSlowLate - groupFastLate,
  SlowLateVsFastEarly = groupSlowLate - groupFastEarly,
  SlowLateVsSlowEarly = groupSlowLate - groupSlowEarly,
  SlowEarlyVsFastLate = groupSlowEarly - groupFastLate,
  SlowEarlyVsSlowLate = groupSlowEarly - groupSlowLate,
  SlowEarlyVsFastEarly = groupSlowEarly - groupFastEarly,
  FastLateVsFastEarly = groupFastLate - groupFastEarly,
  FastLateVsSlowEarly = groupFastLate - groupSlowEarly,
  FastLateVsSlowLate = groupFastLate - groupSlowLate,
  FastEarlyVsFastLate = groupFastEarly - groupFastLate,
  FastEarlyVsSlowLate = groupFastEarly - groupSlowLate,
  FastEarlyVsSlowEarly = groupFastEarly - groupSlowEarly,
  levels = design
)

performDifferentialExpression(nonControlEx,
                              design,
                              cont.matrix,
                              experimentName,
                              nonControlClinical,
                              TRUE)

### Site of Onset Experiment
clinical$group <- clinical$sightOnsetExperiment
nonControlClinical$group <- nonControlClinical$sightOnsetExperiment

experimentName <- "SightOfOnsetVsControl"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + batchControl
  ,
  clinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  ControlVsBulbarLate =   groupControl - groupBulbarLate,
  ControlVsBulbarEarly =   groupControl - groupBulbarEarly,
  ControlVsLimbLate =   groupControl - groupLimbLate,
  ControlVsLimbEarly =   groupControl - groupLimbEarly,
  BulbarLateVsControl =   groupBulbarLate - groupControl,
  BulbarEarlyVsControl =   groupBulbarEarly - groupControl,
  LimbLateVsControl =   groupLimbLate - groupControl,
  LimbEarlyVsControl =   groupLimbEarly - groupControl,
  levels = design
)

performDifferentialExpression(ex, design, cont.matrix, experimentName)

experimentName <- "BulbarVsLimb"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + timeFromOnsetToVisitInYears
  + timeFromVisit1InYears
  + fastSlow
  + batchControl
  ,
  nonControlClinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  BulbarEarlyVsBulbarLate =  groupBulbarEarly - groupBulbarLate,
  BulbarEarlyVsLimbEarly = groupBulbarEarly - groupLimbEarly,
  BulbarEarlyVsLimbLate = groupBulbarEarly - groupLimbLate,
  BulbarLateVsBulbarEarly = groupBulbarLate - groupBulbarEarly,
  BulbarLateVsLimbEarly = groupBulbarLate - groupLimbEarly,
  BulbarLateVsLimbLate = groupBulbarLate - groupLimbLate,
  LimbEarlyVsBulbarEarly = groupLimbEarly - groupBulbarEarly,
  LimbEarlyVsBulbarLate = groupLimbEarly - groupBulbarLate,
  LimbEarlyVsLimbLate = groupLimbEarly - groupLimbLate,
  LimbLateVsBulbarEarly = groupLimbLate - groupBulbarEarly,
  LimbLateVsBulbarLate = groupLimbLate - groupBulbarLate,
  LimbLateVsLimbEarly =   groupLimbLate - groupLimbEarly,
  levels = design
)

performDifferentialExpression(nonControlEx,
                              design,
                              cont.matrix,
                              experimentName,
                              nonControlClinical,
                              TRUE)


### Progression and Site of Onset Experiment
clinical$group <- clinical$progressionAndSightOnsetExperiment
nonControlClinical$group <- nonControlClinical$progressionAndSightOnsetExperiment

experimentName <- "ProgressionAndSightOfOnsetVsControl"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + batchControl
  ,
  clinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  ControlVsFastBulbarEarly =   groupControl - groupFastBulbarEarly,
  ControlVsFastBulbarLate =   groupControl - groupFastBulbarLate,
  ControlVsFastLimbEarly =   groupControl - groupFastLimbEarly,
  ControlVsFastLimbLate =   groupControl - groupFastLimbLate,
  ControlVsSlowBulbarEarly =   groupControl - groupSlowBulbarEarly,
  ControlVsSlowBulbarLate =   groupControl - groupSlowBulbarLate,
  ControlVsSlowLimbEarly =   groupControl - groupSlowLimbEarly,
  ControlVsSlowLimbLate =   groupControl - groupSlowLimbLate,
  FastBulbarEarlyVsControl =   groupFastBulbarEarly - groupControl,
  FastBulbarLateVsControl =   groupFastBulbarLate - groupControl,
  FastLimbEarlyVsControl =   groupFastLimbEarly - groupControl,
  FastLimbLateVsControl =   groupFastLimbLate - groupControl,
  SlowBulbarEarlyVsControl =   groupSlowBulbarEarly - groupControl,
  SlowBulbarLateVsControl =   groupSlowBulbarLate - groupControl,
  SlowLimbEarlyVsControl =   groupSlowLimbEarly - groupControl,
  SlowLimbLateVsControl =   groupSlowLimbLate - groupControl,
  levels = design
)

performDifferentialExpression(ex, design, cont.matrix, experimentName)

experimentName <- "ProgressionAndSightOfOnset"

design <- model.matrix(
  ~ 0
  + group
  + ageAtVisit
  + ethnicity
  + gender
  + sampleStorageDays
  + statinUse
  + timeFromOnsetToVisitInYears
  + timeFromVisit1InYears
  #+ fastSlow
  #+ BulbarLimb
  + batchControl
  ,
  nonControlClinical
)

# set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  groupFastBulbarEarlyVsgroupFastBulbarLate = groupFastBulbarEarly - groupFastBulbarLate,
  groupFastBulbarEarlyVsgroupFastLimbEarly = groupFastBulbarEarly - groupFastLimbEarly,
  groupFastBulbarEarlyVsgroupFastLimbLate = groupFastBulbarEarly - groupFastLimbLate,
  groupFastBulbarEarlyVsgroupSlowBulbarEarly = groupFastBulbarEarly - groupSlowBulbarEarly,
  groupFastBulbarEarlyVsgroupSlowBulbarLate = groupFastBulbarEarly - groupSlowBulbarLate,
  groupFastBulbarEarlyVsgroupSlowLimbEarly = groupFastBulbarEarly - groupSlowLimbEarly,
  groupFastBulbarEarlyVsgroupSlowLimbLate = groupFastBulbarEarly - groupSlowLimbLate,
  groupFastBulbarLateVsgroupFastBulbarEarly = groupFastBulbarLate - groupFastBulbarEarly,
  groupFastBulbarLateVsgroupFastLimbEarly = groupFastBulbarLate - groupFastLimbEarly,
  groupFastBulbarLateVsgroupFastLimbLate = groupFastBulbarLate - groupFastLimbLate,
  groupFastBulbarLateVsgroupSlowBulbarEarly = groupFastBulbarLate - groupSlowBulbarEarly,
  groupFastBulbarLateVsgroupSlowBulbarLate = groupFastBulbarLate - groupSlowBulbarLate,
  groupFastBulbarLateVsgroupSlowLimbEarly = groupFastBulbarLate - groupSlowLimbEarly,
  groupFastBulbarLateVsgroupSlowLimbLate = groupFastBulbarLate - groupSlowLimbLate,
  groupFastLimbEarlyVsgroupFastBulbarEarly = groupFastLimbEarly - groupFastBulbarEarly,
  groupFastLimbEarlyVsgroupFastBulbarLate = groupFastLimbEarly - groupFastBulbarLate,
  groupFastLimbEarlyVsgroupFastLimbLate = groupFastLimbEarly - groupFastLimbLate,
  groupFastLimbEarlyVsgroupSlowBulbarEarly = groupFastLimbEarly - groupSlowBulbarEarly,
  groupFastLimbEarlyVsgroupSlowBulbarLate = groupFastLimbEarly - groupSlowBulbarLate,
  groupFastLimbEarlyVsgroupSlowLimbEarly = groupFastLimbEarly - groupSlowLimbEarly,
  groupFastLimbEarlyVsgroupSlowLimbLate = groupFastLimbEarly - groupSlowLimbLate,
  groupFastLimbLateVsgroupFastBulbarEarly = groupFastLimbLate - groupFastBulbarEarly,
  groupFastLimbLateVsgroupFastBulbarLate = groupFastLimbLate - groupFastBulbarLate,
  groupFastLimbLateVsgroupFastLimbEarly = groupFastLimbLate - groupFastLimbEarly,
  groupFastLimbLateVsgroupSlowBulbarEarly = groupFastLimbLate - groupSlowBulbarEarly,
  groupFastLimbLateVsgroupSlowBulbarLate = groupFastLimbLate - groupSlowBulbarLate,
  groupFastLimbLateVsgroupSlowLimbEarly = groupFastLimbLate - groupSlowLimbEarly,
  groupFastLimbLateVsgroupSlowLimbLate = groupFastLimbLate - groupSlowLimbLate,
  groupSlowBulbarEarlyVsgroupFastBulbarEarly = groupSlowBulbarEarly - groupFastBulbarEarly,
  groupSlowBulbarEarlyVsgroupFastBulbarLate = groupSlowBulbarEarly - groupFastBulbarLate,
  groupSlowBulbarEarlyVsgroupFastLimbEarly = groupSlowBulbarEarly - groupFastLimbEarly,
  groupSlowBulbarEarlyVsgroupFastLimbLate = groupSlowBulbarEarly - groupFastLimbLate,
  groupSlowBulbarEarlyVsgroupSlowBulbarLate = groupSlowBulbarEarly - groupSlowBulbarLate,
  groupSlowBulbarEarlyVsgroupSlowLimbEarly = groupSlowBulbarEarly - groupSlowLimbEarly,
  groupSlowBulbarEarlyVsgroupSlowLimbLate = groupSlowBulbarEarly - groupSlowLimbLate,
  groupSlowBulbarLateVsgroupFastBulbarEarly = groupSlowBulbarLate - groupFastBulbarEarly,
  groupSlowBulbarLateVsgroupFastBulbarLate = groupSlowBulbarLate - groupFastBulbarLate,
  groupSlowBulbarLateVsgroupFastLimbEarly = groupSlowBulbarLate - groupFastLimbEarly,
  groupSlowBulbarLateVsgroupFastLimbLate = groupSlowBulbarLate - groupFastLimbLate,
  groupSlowBulbarLateVsgroupSlowBulbarEarly = groupSlowBulbarLate - groupSlowBulbarEarly,
  groupSlowBulbarLateVsgroupSlowLimbEarly = groupSlowBulbarLate - groupSlowLimbEarly,
  groupSlowBulbarLateVsgroupSlowLimbLate = groupSlowBulbarLate - groupSlowLimbLate,
  groupSlowLimbEarlyVsgroupFastBulbarEarly = groupSlowLimbEarly - groupFastBulbarEarly,
  groupSlowLimbEarlyVsgroupFastBulbarLate = groupSlowLimbEarly - groupFastBulbarLate,
  groupSlowLimbEarlyVsgroupFastLimbEarly = groupSlowLimbEarly - groupFastLimbEarly,
  groupSlowLimbEarlyVsgroupFastLimbLate = groupSlowLimbEarly - groupFastLimbLate,
  groupSlowLimbEarlyVsgroupSlowBulbarEarly = groupSlowLimbEarly - groupSlowBulbarEarly,
  groupSlowLimbEarlyVsgroupSlowBulbarLate = groupSlowLimbEarly - groupSlowBulbarLate,
  groupSlowLimbEarlyVsgroupSlowLimbLate = groupSlowLimbEarly - groupSlowLimbLate,
  groupSlowLimbLateVsgroupFastBulbarEarly = groupSlowLimbLate - groupFastBulbarEarly,
  groupSlowLimbLateVsgroupFastBulbarLate = groupSlowLimbLate - groupFastBulbarLate,
  groupSlowLimbLateVsgroupFastLimbEarly = groupSlowLimbLate - groupFastLimbEarly,
  groupSlowLimbLateVsgroupFastLimbLate = groupSlowLimbLate - groupFastLimbLate,
  groupSlowLimbLateVsgroupSlowBulbarEarly = groupSlowLimbLate - groupSlowBulbarEarly,
  groupSlowLimbLateVsgroupSlowBulbarLate = groupSlowLimbLate - groupSlowBulbarLate,
  groupSlowLimbLateVsgroupSlowLimbEarly = groupSlowLimbLate - groupSlowLimbEarly,
  levels = design
)

performDifferentialExpression(nonControlEx,
                              design,
                              cont.matrix,
                              experimentName,
                              nonControlClinical,
                              TRUE)
