library(limma)
ex <-
  read.csv("data/lipidomics/expressionDataRaw.csv", row.names = 1)

ex <- ex[names(sort(sapply(ex, max, na.rm = T), decreasing = T))]

asinhTransformedEx <- asinh(ex)

write.csv(asinhTransformedEx,
          "data/lipidomics/asinhTransformedExpressionDataAllSamples.csv")

outliers <- c("#FE7", #"FL7",
              "FL10", "FE10")

ex <- ex[,-which(names(ex) %in% outliers)]

ex <- asinh(ex)

write.csv(
  ex,
  "data/lipidomics/asinhTransformedExpressionDataOutliersAndDuplicatesRemoved.csv"
)

outliers <-
  append(
    outliers,
    c(
      "H9",
      "H8",
      "H5",
      "H2"
    )
  )

ex <- ex[,-which(names(ex) %in% outliers)]

ex <- normalizeBetweenArrays(ex) # normalize data

write.csv(
  ex,
  "data/lipidomics/normalisedasinhTransformedExpressionDataRawOutliersAndDuplicatesRemoved.csv"
)


#####
ex <-
  read.csv("data/lipidomics/expressionDataRawIncludingTotals.csv", row.names = 1)

ex <- ex[names(sort(sapply(ex, max, na.rm = T), decreasing = T))]

asinhTransformedEx <- asinh(ex)

write.csv(asinhTransformedEx,
          "data/lipidomics/asinhTransformedexpressionDataRawIncludingTotals.csv")

outliers <- c("#FE7", #"FL7",
              "FL10", "FE10")

ex <- ex[,-which(names(ex) %in% outliers)]

ex <- asinh(ex)

write.csv(
  ex,
  "data/lipidomics/asinhTransformedexpressionDataRawIncludingTotalsOutliersAndDuplicatesRemoved.csv"
)

outliers <-
  append(
    outliers,
    c(
      "H9",
      "H8",
      "H5",
      "H2"
    )
  )

ex <- ex[,-which(names(ex) %in% outliers)]

ex <- normalizeBetweenArrays(ex) # normalize data

write.csv(
  ex,
  "data/lipidomics/normalisedasinhTransformedexpressionDataRawIncludingTotalsRawOutliersAndDuplicatesRemoved.csv"
)

