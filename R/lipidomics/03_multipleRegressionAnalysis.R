library("heatmaply")

ex <- read.csv(
  "data/lipidomics/normalisedasinhTransformedExpressionDataRawOutliersAndDuplicatesRemoved.csv",
  row.names = 1
)

ex <- na.omit(ex)

transposedEx <- t(ex)

clinical <-
  read.csv("data/lipidomics/clinicalData.csv", row.names = 1)

head(clinical)

combinedDataset <-
  merge(transposedEx, clinical, by = 0, all.x = TRUE)

combinedDataset <-
  combinedDataset[!is.na(combinedDataset$survivalInMonthsFromVisitDate), ]

colnames(combinedDataset) <-
  gsub("-", "", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub(" ", "", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub(",", "", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[0]", "Zero", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[1]", "One", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[2]", "Two", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[3]", "Three", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[4]", "Four", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[5]", "Five", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[6]", "Six", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[7]", "Seven", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[8]", "Eight", colnames(combinedDataset))
colnames(combinedDataset) <-
  gsub("[9]", "Nine", colnames(combinedDataset))

correlationsToSurvival <- cor(x = combinedDataset[,colnames(combinedDataset)[2:nrow(ex)]],
              y = combinedDataset[,"survivalInMonthsFromVisitDate"])

correlationsToSurvival[order(correlationsToSurvival),]

### Linear model
formulaString <-
  paste("survivalInMonthsFromVisitDate ~ ",
        paste(colnames(combinedDataset)[2:10
                                        #nrow(ex)
                                        ]
                                        , collapse = ' + '),
        collapse = "")
formulaString <-
  as.formula(formulaString)

model <-
  lm(formulaString, data = combinedDataset)
summary(model)


# Mixed Linear model
library(nlme)
fit1 <-
  lme(fixed = formulaString,
      random = ~ 1 | Patient.ID,
      data = combinedDataset)
summary(fit1)
