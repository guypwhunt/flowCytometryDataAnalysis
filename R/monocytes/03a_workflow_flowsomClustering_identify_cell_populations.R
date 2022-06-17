try(source("R/01_functions.R"))

loadlibraries()

directory <- "monocytes"

queries <- list(
  c("CD14...BV605.A" = "low",
    "CD16...PE.CF595.A" = "high",
    "CD11b...17BV421.A" = "low",
    "CD11b.activated...PE.Cy7.A" = "low"),
  c("CD14...BV605.A" = "low",
    "CD16...PE.CF595.A" = "high",
    "CD11b...17BV421.A" = "high"),
  c("CD14...BV605.A" = "low",
    "CD16...PE.CF595.A" = "high",
    "CD11b.activated...PE.Cy7.A" = "high"),

  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "high",
    "CD11b...17BV421.A" = "low",
    "CD11b.activated...PE.Cy7.A" = "low"),
  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "high",
    "CD11b...17BV421.A" = "high"),
  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "high",
    "CD11b.activated...PE.Cy7.A" = "high"),

  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "low",
    "CD11b...17BV421.A" = "low",
    "CD11b.activated...PE.Cy7.A" = "low"),
  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "low",
    "CD11b...17BV421.A" = "high"),
  c("CD14...BV605.A" = "high",
    "CD16...PE.CF595.A" = "low",
    "CD11b.activated...PE.Cy7.A" = "high"))

names(queries) <- c("Non-Classical Monocytes",
                    "CD11b+ Non-Classical Monocytes",
                    "CD11b+ Non-Classical Monocytes",

                    "Intermediate Monocytes",
                    "CD11b+ Intermediate Monocytes",
                    "CD11b+ Intermediate Monocytes",

                    "Classical Monocytes",
                    "CD11b+ Classical Monocytes",
                    "CD11b+ Classical Monocytes")

defineFlowSomCellPopulations(directory, queries)
