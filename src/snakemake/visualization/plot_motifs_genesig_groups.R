invisible(sapply(c(
  "stringr",
  "progress",
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

sobj@meta.data$genesig_group <- "Q"
sobj@meta.data[
  sobj[[]]$is_earlydiff == "earlydiff" &
    sobj[[]]$is_HSC == "HSC",
][["genesig_group"]] <- "B"
sobj@meta.data[
  sobj[[]]$is_earlydiff == "not_earlydiff" &
    sobj[[]]$is_HSC == "HSC",
][["genesig_group"]] <- "H"
sobj@meta.data[
  sobj[[]]$is_earlydiff == "earlydiff" &
    sobj[[]]$is_HSC == "not_HSC",
][["genesig_group"]] <- "E"
sobj@meta.data[
  sobj[[]]$is_earlydiff == "not_earlydiff" &
    sobj[[]]$is_HSC == "not_HSC",
][["genesig_group"]] <- "N"

# table(sobj[[c("is_earlydiff", "is_HSC", "genesig_group")]])

################################################################################
#                              Plotting functions                              #
################################################################################

##################################################
#  DimPlots to show location of the four groups  #
##################################################

plot_genesig_group <- function(genesig_group, umap = "ATAC") {
  Idents(sobj) <- "genesig_group"
  if (genesig_group == "N") {
    group <- "HSC: [ ]\nEarlydiff: [ ]"
  } else if (genesig_group == "H") {
    group <- "HSC: [x]\nEarlydiff: [ ]"
  } else if (genesig_group == "E") {
    group <- "HSC: [ ]\nEarlydiff: [x]"
  } else if (genesig_group == "B") {
    group <- "HSC: [x]\nEarlydiff: [x]"
  }
  axis_text_style <- element_text(
    face = "italic",
    size = 10
  )
  p <- DimPlot(sobj,
    cells.highlight = WhichCells(sobj, idents = genesig_group),
    reduction = str_c("UMAPharmony", umap)
  ) +
    coord_fixed() +
    labs(
      title = group,
      subtitle = str_c(
        length(sobj[["genesig_group"]][
          sobj[["genesig_group"]]$genesig_group ==
            genesig_group,
        ]),
        " cells"
      )
    ) +
    NoLegend() +
    theme(
      plot.title = element_text(
        face = "plain",
        size = 10
      ),
      plot.subtitle = element_text(
        face = "plain",
        size = 8
      ),
      axis.title = element_blank(),
      axis.text.x = axis_text_style,
      axis.text.y = axis_text_style
    )
  return(p)
}

compile_genesig_group <- function(umap = "ATAC") {
  p <- (
    plot_genesig_group("B", umap = umap) |
      plot_genesig_group("H", umap = umap)
  ) / (
    plot_genesig_group("E", umap = umap) |
      plot_genesig_group("N", umap = umap)
  ) +
    plot_annotation(
      title = "Gene signature-classified groups on UMAP",
      caption = str_c("Harmony-integrated UMAP for ", umap, " data"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          hjust = .5
        )
      )
    ) &
    theme(
      plot.caption = element_text(
        face = "italic",
        size = 8
      )
    )
  return(p)
}

###################################################################
#  FeaturePlots to show motif activity score per group and assay  #
###################################################################

plot_genesig_motif <- function(assay,
                               motif_idx,
                               genesig_group,
                               umap = "ATAC") {
  Idents(sobj) <- "genesig_group"
  if (genesig_group == "N") {
    group <- "HSC: [ ]\nEarlydiff: [ ]"
  } else if (genesig_group == "H") {
    group <- "HSC: [x]\nEarlydiff: [ ]"
  } else if (genesig_group == "E") {
    group <- "HSC: [ ]\nEarlydiff: [x]"
  } else if (genesig_group == "B") {
    group <- "HSC: [x]\nEarlydiff: [x]"
  }
  axis_text_style <- element_text(
    face = "italic",
    size = 10
  )
  p <- FeaturePlot(sobj,
    features = str_c(assay, "chromvar_", motif_ids[motif_idx]),
    reduction = str_c("UMAPharmony", umap),
    cells = WhichCells(sobj, idents = genesig_group),
    min.cutoff = "q1", max.cutoff = "q99",
    pt.size = .5,
    order = T
  ) +
    theme(axis.title.y.right = element_blank()) +
    coord_fixed() +
    NoLegend() +
    NoAxes() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    )
  if (assay == "proximal") {
    p <- p +
      labs(
        title = group,
        subtitle = str_c(
          length(sobj[["genesig_group"]][
            sobj[["genesig_group"]]$genesig_group ==
              genesig_group,
          ]),
          " cells"
        )
      ) +
      theme(
        plot.title = element_text(
          face = "plain",
          hjust = .5,
          size = 10
        ),
        plot.subtitle = element_text(
          face = "plain",
          hjust = .5,
          size = 8
        )
      )
  } else if (assay == "distal") {
    p <- p +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank()
      )
  }
  return(p)
}
# plot_genesig_motif(assay = "proximal", motif_idx = 4, genesig_group = "B")

plot_assay_genesig_motif <- function(motif_idx,
                                     assay = "proximal",
                                     umap = "ATAC") {
  p <-
    plot_genesig_motif(
      genesig_group = "B",
      assay = assay,
      motif_idx = motif_idx,
      umap = umap
    ) |
      plot_genesig_motif(
        genesig_group = "H",
        assay = assay,
        motif_idx = motif_idx,
        umap = umap
      ) |
      plot_genesig_motif(
        genesig_group = "E",
        assay = assay,
        motif_idx = motif_idx,
        umap = umap
      ) |
      (
        plot_genesig_motif(
          genesig_group = "N",
          assay = assay,
          motif_idx = motif_idx,
          umap = umap
        ) +
          RestoreLegend() +
          labs(
            color = assay
          ) +
          theme(
            legend.key.size = unit(.5, "cm"),
            legend.title = element_text(
              size = 10
            ),
            legend.text = element_text(
              size = 8
            )
          )
      )
  p <- p +
    theme(
      axis.title.y.right = element_text(
        face = "plain"
      )
    )
  return(p)
}
# plot_assay_genesig_motif(4)

compile_genesig_motif <- function(motif_idx,
                                  umap = "ATAC") {
  p <- (
    plot_assay_genesig_motif(
      assay = "proximal",
      motif_idx = motif_idx,
      umap = umap
    )
  ) / (
    plot_assay_genesig_motif(
      assay = "distal",
      motif_idx = motif_idx,
      umap = umap
    )
  ) +
    plot_annotation(
      title = str_c(
        "Cells by gene signature-classified groups on UMAP,\n",
        "coloured by ",
        sobj[["proximal"]]@motifs@motif.names[motif_ids[motif_idx]],
        " (",
        motif_ids[motif_idx],
        ") motif activity score"
      ),
      caption = str_c("Harmony-integrated UMAP for ", umap, " data"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          hjust = .5
        )
      )
    ) &
    theme(
      plot.caption = element_text(
        face = "italic",
        size = 8
      )
    )
  return(p)
}
# compile_genesig_motif(4)

################################################################################
#                     Carry out plotting and save results                      #
################################################################################

for (modality in c("ATAC", "RNA")) {
  ggsave(
    compile_genesig_group(umap = modality),
    filename = snakemake@output[[str_c(
      "plot_umap_",
      str_to_lower(modality)
    )]],
    device = "svg",
    units = "in",
    width = 10,
    height = 10,
    dpi = "retina"
  )
}

# compile_genesig_group(umap = "ATAC")

n_iter <- length(motif_ids)
pb <- progress_bar$new(
  format = str_c(
    "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || ",
    "Estimated time remaining: :eta]"
  ),
  total = n_iter * 2,
  complete = "=", # Completion bar character
  incomplete = "-", # Incomplete bar character
  current = ">", # Current bar character
  clear = FALSE, # If TRUE, clears the bar when finish
  width = 100 # Width of the progress bar
)
for (idx in seq_len(n_iter)) {
  for (modality in c("ATAC", "RNA")) {
    p <- compile_genesig_motif(idx, umap = modality)
    ggsave(
      plot = p,
      filename = str_c(
        snakemake@output[["plot_dir"]], "/",
        motif_ids[idx], "_",
        motif_names[idx],
        "_",
        modality,
        ".svg"
      ),
      units = "in",
      width = 11,
      height = 7,
      dpi = "retina",
      device = "svg"
    )
    pb$tick()
  }
}

# # lapply(seq_len(length(motif_ids)), FUN = function(idx) {
# sink <- lapply(seq_len(5), FUN = function(idx) {
#   for (modality in c("ATAC", "RNA")) {
#     p <- compile_genesig_motif(idx, umap = modality)
#     ggsave(
#       plot = p,
#       filename = str_c(
#         snakemake@output[["plot_dir"]], "/",
#         motif_ids[idx], "_",
#         motif_names[idx],
#         "_",
#         modality,
#         ".svg"
#       ),
#       units = "in",
#       width = 10,
#       height = 6,
#       dpi = "retina",
#       device = "svg"
#     )
#   }
# })
