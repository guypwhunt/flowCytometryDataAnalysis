try(source("R/01_functions.R"))
try(source("R/00_datasets.R"))

loadlibraries()

directoryNames <- c(
  # "gpr18BCells",
  # "gpr18Monocytes",
  # "gpr18Senescence",
  # "gpr18TCells",
  # "gpr32BCells"#,
  "gpr32Monocytes",
  # "gpr32Senescence",
  # "gpr32TCells"
)

clusterNames <- clusterColumns[4]

# clusterName <- clusterNames[2]
#
# directoryName <- "gpr32BCells"

for (directoryName in directoryNames) {
  message()
  message(directoryName)
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

      df <- df[order(df$typeOfCells),]

      df$ID <- as.factor(as.numeric(factor(df$typeOfCells, levels = unique(df$typeOfCells))))

      df$typeOfCells <-
        paste0(df$typeOfCells, " (", df$ID, ")")

      df$typeOfCells <-
        factor(df$typeOfCells, levels = unique(df$typeOfCells))


      label <-
        df %>% group_by(clusters_phenograph, typeOfCells) %>%
        dplyr::select(umap_1, umap_2) %>% summarize_all(mean) %>%
        as.data.frame()

      par(mar = c(1, 1, 1, 1))

      plot <-
        ggplot(df, aes(
          x = umap_1,
          y = umap_2,
          color = as.factor(clusters_phenograph)
        )) +
        geom_point(size = 0.1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_label_repel(
          aes(label = as.integer(clusters_phenograph)),
          data = label
          ,
          size = 2,
          force_pull = 0,
          max.time = 2,
          show.legend = FALSE
        ) +
        labs(y = "UMAP 2", x = "UMAP 1")  +
        guides(
          color = guide_legend(
            title = "Cell Populations",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1)          ),
          fill = guide_legend(
            title = "Cell Populations",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1)          ),
          guide = guide_legend(
            title = "Clusters",
            override.aes = list(shape = 16, size = 3,
                                ncol = 1)
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

