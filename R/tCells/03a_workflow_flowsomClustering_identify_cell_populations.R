try(source("R/01_functions.R"))

loadlibraries()

directory <- "senescence"

queries <- list(
  "Viral Associated Senescent CD8+ T Cells" = c(
    "CD8.PerCP.Cy5.5.A" = "high",
    "CD27.BV421.A" = "low",
    "CD45RA.BV605.A" = "high",
    "KLRG1.PE.A" = "high"
  ),
  "Non-Viral Associated Senescent CD8+ T Cells" = c(
    "CD8.PerCP.Cy5.5.A" = "high",
    "CD27.BV421.A" = "low",
    "CD45RA.BV605.A" = "high",
    "KLRG1.PE.A" = "low"
  ),

  "Intermediate Senescent 2 CD8+ T Cells" =
    c(
      "CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "high"
    ),
  "Late Senescent CD8+ T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "low"),
  "Intermediate Senescent 1 CD8+ T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "low"),
  "Early Senescent CD8+ T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "high"),

  "Viral Associated Senescent CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CD27.BV421.A" = "low",
      "CD45RA.BV605.A" = "high",
      "KLRG1.PE.A" = "high"),
  "Non-Viral Associated Senescent CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CD27.BV421.A" = "low",
      "CD45RA.BV605.A" = "high",
      "KLRG1.PE.A" = "low"),

  "Intermediate Senescent 2 CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "high"),
  "Late Senescent CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "low"),
  "Intermediate Senescent 1 CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "low"),
  "Early Senescent CD4+ T Cells" =
    c("CD4.PE.CF594.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "high"),

  "Naive CD8+ T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CD45RA.BV605.A" = "low"),
  "Naive CD4+ T Cells",
  c("CD45RA.BV605.A" = "low",
    "CD4.PE.CF594.A" = "low"),

  "Double-Postive T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CD4.PE.CF594.A" = "high"),
  "Double-Negative T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "low",
      "CD4.PE.CF594.A" = "low"),

  "Naive Double-Postive T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "high",
      "CD4.PE.CF594.A" = "high",
      "CD45RA.BV605.A" = "low"),

  "Naive Double-Negative T Cells" =
    c("CD8.PerCP.Cy5.5.A" = "low",
      "CD4.PE.CF594.A" = "low",
      "CD45RA.BV605.A" = "low"),

  "Viral Associated Senescent Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CD27.BV421.A" = "low",
      "CD45RA.BV605.A" = "high",
      "KLRG1.PE.A" = "high"),

  "Non-Viral Associated Senescent Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CD27.BV421.A" = "low",
      "CD45RA.BV605.A" = "high",
      "KLRG1.PE.A" = "low"),

  "Intermediate Senescent 2 Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "high"),
  "Late Senescent Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "low",
      "CD28.BV785.A" = "low"),
  "Intermediate Senescent 1 Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "low"),
  "Early Senescent Double-Postive T Cells" =
    c("CD4.PE.CF594.A" = "high", "CD8.PerCP.Cy5.5.A" = "high",
      "CCR7.PE.Cy7.A" = "low",
      "CD45RA.BV605.A" = "high",
      "CD27.BV421.A" = "high",
      "CD28.BV785.A" = "high")
)

defineFlowSomCellPopulations(directory, queries)

