invisible(sapply(c(
  "stringr",
  "patchwork",
  "ggplot2",
  "Seurat",
  "GenomicRanges",
  "Signac"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

dir.create(snakemake@output[["plot_dir"]], recursive = TRUE)
stopifnot(file.exists(snakemake@output[["plot_dir"]]))

################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])
DefaultAssay(sobj) <- "proximal"
Idents(sobj) <- "origin"

motif_names <- ConvertMotifID(sobj, id = colnames(Motifs(sobj[["proximal"]])))
motif_ids <- colnames(Motifs(sobj[["distal"]]))

################################################################################
#                              Plotting functions                              #
################################################################################

plot_distal_proximal_split <- function(assay, motif_idx, split = "origin") {
  save_assay <- DefaultAssay(sobj)
  Seurat::DefaultAssay(sobj) <- assay
  plot_title <- str_c(
    "Motif activity score for ", assay, " peaks"
  )
  plot_subtitle <- str_c(
    "ID: \t\t\t", motif_ids[motif_idx], "\n",
    "Name: \t", sobj[[assay]]@motifs@motif.names[motif_ids[motif_idx]]
  )
  axis_text_style <- element_text(
    face = "italic",
    size = 10
  )
  p <- FeaturePlot(sobj,
    features = str_c(assay, "chromvar_", motif_ids[motif_idx]),
    reduction = "UMAPharmonyATAC",
    min.cutoff = "q1", max.cutoff = "q99",
    order = T,
    split.by = split
    # combine = F
  )
  # ) +
  #   labs(
  #     title = plot_title,
  #     subtitle = plot_subtitle,
  #     # tag = str_c("[", str_sub(assay, end = 4), "]"),
  #     caption = "Harmony-integrated UMAP for ATAC data"
  #   ) +
  #   xlab("") +
  #   ylab("") +
  #   theme(
  #     axis.text.x = axis_text_style,
  #     axis.text.y = axis_text_style,
  #     plot.caption = element_text(
  #       face = "italic"
  #     ),
  #     plot.tag = element_text(
  #       face = "plain",
  #       size = "12"
  #     )
  # )
  # p[[1]] <- p[[1]] + coord_fixed()
  # p[[2]] <- p[[2]] + coord_fixed()
  p <- p + theme(axis.title.y.right = element_blank())
  p <- lapply(list(p[[1]], p[[2]]), FUN = function(plot) {
    return(
      plot +
        coord_fixed() +
        xlab("") +
        ylab("") +
        theme(
          axis.text.x = axis_text_style,
          axis.text.y = axis_text_style,
          plot.caption = element_text(
            face = "italic"
          ),
          plot.tag = element_text(
            face = "plain",
            size = "12"
          )
        )
    )
  })
  p <- p[[1]] | (p[[2]] + RestoreLegend())
  Seurat::DefaultAssay(sobj) <- save_assay
  return(p)
}

compile_final_plot <- function(motif_idx) {
  pprox <- plot_distal_proximal_split("proximal", motif_idx)
  pdist <- plot_distal_proximal_split("distal", motif_idx)
  plot_subtitle <- str_c(
    "ID: \t\t\t", motif_ids[motif_idx], "\n",
    "Name: \t", sobj[["proximal"]]@motifs@motif.names[motif_ids[motif_idx]]
  )
  final_p <- (
    pprox[[1]] +
      ylab("proximal peaks") +
      theme(axis.title = element_text(
        face = "bold",
        size = 10
      ))
  ) +
    pprox[[2]] +
    (
      pdist[[1]] +
        ggtitle("") +
        ylab("distal peaks") +
        theme(axis.title = element_text(
          face = "bold",
          size = 10
        ))
    ) +
    (
      pdist[[2]] + ggtitle("")
    ) +
    plot_layout(nrow = 2, ncol = 2) +
    plot_annotation(
      title = "Motif activity scores, split by sample",
      subtitle = plot_subtitle,
      caption = "Harmony-integrated UMAP for ATAC data"
    ) &
    theme(
      title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        face = "plain",
        size = 12
      ),
      plot.caption = element_text(
        face = "italic",
        size = 10
      )
    )
  return(final_p)
}

################################################################################
#                     Carry out plotting and save results                      #
################################################################################

lapply(seq_len(length(motif_ids)), FUN = function(idx) {
# lapply(seq_len(10), FUN = function(idx) {
  p <- compile_final_plot(idx)
  ggsave(
    plot = p,
    # filename = str_c("data/adhoc/ATAC/signac-motifs/motif_", idx, ".svg"),
    filename = str_c(
      snakemake@output[["plot_dir"]], "/",
      motif_ids[idx], "_",
      motif_names[idx],
      ".svg"
    ),
    units = "in",
    width = 10,
    height = 10,
    dpi = "retina",
    device = "svg"
  )
})
