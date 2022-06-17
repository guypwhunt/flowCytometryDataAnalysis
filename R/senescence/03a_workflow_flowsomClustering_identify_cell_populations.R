try(source("R/01_functions.R"))

loadlibraries()

directory <- "tCells"

queries <- list(
  c("CD8.BV650.A" = "high",
    "CD45RO.PE.CF595.A" = "high"),
  c("CD8.BV650.A" = "high",
    "CD45RO.PE.CF595.A" = "low"),

  c("CD4.PerCP.Cy5.5.A" = "high",
    "CD45RO.PE.CF595.A" = "high"),
  c("CD4.PerCP.Cy5.5.A" = "high",
    "CD45RO.PE.CF595.A" = "low"),

  c("CD4.PerCP.Cy5.5.A" = "high",
    "CD45RO.PE.CF595.A" = "high",
    "CD25.BV786.A" = "high",
    "CD127.BV510.A" = "low",
    "FoxP3.PE.A" = "high"
  ),
  c("CD4.PerCP.Cy5.5.A" = "high",
    "CD45RO.PE.CF595.A" = "low",
    "CD25.BV786.A" = "high",
    "CD127.BV510.A" = "low",
    "FoxP3.PE.A" = "high"
  ),

  c("CD8.BV650.A" = "high",
    "CD4.PerCP.Cy5.5.A" = "high"),
  c("CD8.BV650.A" = "low",
    "CD4.PerCP.Cy5.5.A" = "low"))

names(queries) <- c("Memory CD8+ T Cells",
                    "Naive CD8+ T Cells",

                    "Memory CD4+ T Cells",
                    "Naive CD4+ T Cells",

                    "Memory T Regulatory Cells",
                    "Naive T Regulatory Cells",

                    "Double-Postive T Cells",
                    "Double-Negative T Cells"
)

defineFlowSomCellPopulations(directory, queries)
