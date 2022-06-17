try(source("R/01_functions.R"))

loadlibraries()

directory <- "bCells"

queries <- list(
  c("CD27...BV650.A" = "high",
    "CD24...BV605.A" = "high",
    "IgD...PerCP.Cy5.5.A" = "high"),
  c("CD27...BV650.A" = "high",
    "CD24...BV605.A" = "low",
    "IgD...PerCP.Cy5.5.A" = "high"),
  c("CD27...BV650.A" = "high",
    "CD24...BV605.A" = "high",
    "IgD...PerCP.Cy5.5.A" = "low"),
  c("CD27...BV650.A" = "high",
    "CD24...BV605.A" = "low",
    "IgD...PerCP.Cy5.5.A" = "low"),
  c("CD27...BV650.A" = "low",
    "CD24...BV605.A" = "low",
    "IgD...PerCP.Cy5.5.A" = "low"),
  c("CD27...BV650.A" = "low",
    "CD24...BV605.A" = "low",
    "IgD...PerCP.Cy5.5.A" = "high"),
  c("CD27...BV650.A" = "low",
    "CD24...BV605.A" = "high",
    "IgD...PerCP.Cy5.5.A" = "low"),
  c("CD27...BV650.A" = "low",
    "CD24...BV605.A" = "high",
    "IgD...PerCP.Cy5.5.A" = "high"))

names(queries) <- c("Unswitched Memory B Cells",
                    "Unswitched Memory B Cells",
                    "Switched Memory B Cells",
                    "Switched Memory B Cells",
                    "Late Memory B Cells",
                    "Naive B Cells",
                    "Immature B Cells",
                    "Follicular B Cells"
)

defineFlowSomCellPopulations(directory, queries)
