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

plot_distal_proximal <- function(assay, motif_idx) {
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
    split.by = NULL
  ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      # tag = str_c("[", str_sub(assay, end = 4), "]"),
      caption = "Harmony-integrated UMAP for ATAC data"
    ) +
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
    ) +
    coord_fixed()
  Seurat::DefaultAssay(sobj) <- save_assay
  return(p)
}
# plot_distal_proximal("proximal", 2) | plot_distal_proximal("distal", 2)

# plot_distal_proximal("proximal", 6) | plot_distal_proximal("distal", 6)

lapply(seq_len(length(motif_ids)), FUN = function(idx) {
# lapply(seq_len(10), FUN = function(idx) {
  p <- (plot_distal_proximal("proximal", idx) |
    plot_distal_proximal("distal", idx))
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
    width = 11,
    height = 6,
    dpi = "retina",
    device = "svg"
  )
})
