try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryNames <- c(
  # "gpr18BCells"#,
  # "gpr18Monocytes"#,
  # "gpr18Senescence",
  "gpr18TCells"#,
  # "gpr32BCells"#,
  # "gpr32Monocytes"#,
  # "gpr32Senescence"#,
  # "gpr32TCells"
)

clusterNames <- clusterColumns[4]

# clusterName <- clusterNames[1]
#
# directoryName <- "gpr32TCells"

for (directoryName in directoryNames) {
  message()
  message(directoryName)

  cellPopulationOrder <- identifyCellPopulationOrder(directoryName)

  for (clusterName in clusterNames) {
    message(clusterName)
    try({
      df <-
        as.data.frame(fread(
          paste0('data/', directoryName, '/clusteringOutput/umapDf.csv')
        ))

      markerPopulations <-
        as.data.frame(fread(
          paste0(
            'data/',
            directoryName,
            '/clusteringOutput/',
            clusterName,
            'Markers.csv'
          )
        ))

      cellPopulations <-
        as.data.frame(fread(
          paste0(
            'data/',
            directoryName,
            '/clusteringOutput/',
            clusterName,
            'CellPopulations.csv'
          )
        ))

      colnames(cellPopulations)[colnames(cellPopulations) == "cell_population"] <-
        "typeOfCells"

      df <-
        merge(
          x = df,
          y = markerPopulations[, c(clusterName, "cell_population")],
          by.x = clusterName,
          by.y = clusterName,
          all.x = TRUE
        )

      df <-
        merge(
          x = df,
          y = cellPopulations[, c(clusterName, "typeOfCells")],
          by.x = clusterName,
          by.y = clusterName,
          all.x = TRUE
        )

      df$cluster_id <- df$cell_population

      df <- updateMarkerNames(df)

      df <- left_join(data.frame(typeOfCells = cellPopulationOrder), df, by = "typeOfCells")

      df <- na.omit(df)

      df$ID <- as.factor(as.numeric(factor(df$typeOfCells, levels = unique(df$typeOfCells))))

      df$typeOfCells <-
        paste0(df$typeOfCells, " (", df$ID, ")")

      df$typeOfCells <-
        factor(df$typeOfCells, levels = unique(df$typeOfCells))

      metacluster_colours <-
        as.vector(colortools::wheel('#F8766D', num = length(levels(df$typeOfCells))))

      label <-
        df %>% group_by(ID, typeOfCells) %>%
        dplyr::select(umap_1, umap_2) %>% summarize_all(mean) %>%
        as.data.frame()

      par(mar = c(1, 1, 1, 1))

      plot <-
        ggplot(df, aes(
          x = umap_1,
          y = umap_2,
          color = as.factor(typeOfCells)
        )) +
        geom_point(size = 0.1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_label_repel(
          aes(label = as.integer(ID)),
          data = label,
          segment.colour = "black",
          size = 4,
          force_pull = 0,
          max.time = 2,
          show.legend = FALSE,
          max.overlaps = Inf,
          fontface = "bold"
        ) +
        labs(y = "UMAP 2", x = "UMAP 1")  +
        scale_color_manual(values = metacluster_colours) +
        guides(
          color = guide_legend(
            title = "Cell Populations",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1),
            ncol = 1
          ),
          fill = guide_legend(
            title = "Cell Populations",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1),
            ncol = 1
          ),
          guide = guide_legend(
            title = "Cell Populations",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1),
            ncol = 1
          )
        )

      print(plot)

    })
  }

  rm(list = ls()[!ls() %in% c("directoryNames", "clusterNames",
                              "directoryName", "clusterName")])
  try(source("R/01_functions.R"))
  try(source("R/00_datasets.R"))
}

