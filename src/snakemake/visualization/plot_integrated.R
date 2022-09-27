suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
set.seed(snakemake@config[["seed"]])

################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])

################################################################################
#                                   Plotting                                   #
################################################################################

single_width <- 16
single_height <- 14
double_width <- 32
double_height <- 16

evil_plot <- function(
    modality = "RNA",
    integration_state = "integrated",
    plot_colouring = "origin"
  ) {
  if (integration_state == "integrated") {
    reduction <- "harmony"
  } else {
    reduction <- "unintegrated"
  }

  if (plot_colouring == "clusters_split") {
    group <- str_c(
      "clusters_",
      modality,
      "_",
      reduction
    )
    p <- DimPlot(
      sobj,
      reduction = str_c("UMAP", reduction, modality),
      group.by = group,
      split.by = "origin",
      shuffle = T
    )

    width <- double_width
    height <- double_height
  } else {
    if (plot_colouring == "clusters") {
      group <- str_c(
        "clusters_",
        modality,
        "_",
        reduction
      )
    } else {
      group <- "origin"
    }

    p <- DimPlot(
      sobj,
      reduction = str_c("UMAP", reduction, modality),
      group.by = group,
      shuffle = T
    )

    width <- single_width
    height <- single_height
  }

  p <- p + ggtitle(str_c(
    str_to_title(integration_state),
    " ",
    modality,
    " coloured on ",
    plot_colouring
  ))

  smk_output <- str_c(
                  "plot_",
                  integration_state,
                  "_",
                  modality,
                  "_",
                  plot_colouring
                )

  ggsave(
    snakemake@output[[smk_output]],
    plot = p,
    device = "svg",
    width = width,
    height = height,
    units = "cm",
    dpi = "retina"
  )
}

modalities <- c("RNA", "ATAC")
integration_states <- c("integrated", "unintegrated")
plot_colourings <- c("origin", "clusters", "clusters_split")

for (modality in modalities) {
  for (integration_state in integration_states) {
    for (plot_colouring in plot_colourings) {
      print(paste0(
        "> Plotting ",
        integration_state,
        " ",
        modality,
        " coloured on ",
        plot_colouring,
        ".."
      ))
      evil_plot(
        modality = modality,
        integration_state = integration_state,
        plot_colouring = plot_colouring
      )
    }
  }
}
