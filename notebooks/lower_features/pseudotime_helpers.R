# invisible(sapply(c(
#   "stringr",
#   "parallel",
#   "docstring",
#   "patchwork",
#   "ggplot2",
#   "tidyr",
#   "dplyr",
#   "tictoc",
#   "readr",
#   "assertthat",
#   "tibble",
#   "slingshot",
#   "ggbeeswarm",
#   "circlize",
#   "ComplexHeatmap",
#   "RColorBrewer",
#   "viridis",
#   "Seurat",
#   "Signac"
# ), FUN = function(x) {
#   suppressPackageStartupMessages(library(x, character.only = T))
# }))

# project_path <- Sys.getenv("PROJECT_PATH")
# output_dir <- str_c(project_path, "/data/adhoc/SAILERX/global_notebooks/")
# plot_dir <- str_c(project_path, "/data/adhoc/SAILERX/global_notebooks/plots/")
# rds_path <- str_c(output_dir, "seurat_object_post_atac_processing.rds")

# set.seed(12345)

# if (!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }
# if (!dir.exists(plot_dir)) {
#   dir.create(plot_dir)
# }

# sailerx_colours <- palette.colors(n = 10, palette = "Paired")
# names(sailerx_colours) <- as.character(0:9)

slingplot <- function(sce,
                      reduction = "UMAP_SAILERX",
                      palette = "viridis",
                      line_type = NULL,
                      subset = FALSE,
                      start_cluster = start_cluster,
                      end_cluster = -1,
                      included_clusters = c(-1),
                      title = "Pseudotime + Trajectories") {
  if (str_to_title(palette) %in% hcl.pals()) {
    colors <- hcl.colors(n = 50, palette = palette)
  }
  if (subset) {
    sce_base <- sce_all
  } else {
    sce_base <- sce
  }
  end_cluster <- end_cluster %>% as.character()
  plot(reducedDims(sce_base)[[reduction]],
    pch = 16, col = rgb(
      red = .6, green = .6, blue = .6, alpha = .7
    ), asp = 1
  )
  points(
    reducedDims(sce)[[reduction]],
    col = colors[
      cut(
        sce[[
          str_c("slingPseudotime_", traj_dict[[end_cluster]][["lvl2"]])
        ]],
        breaks = 50
      )
    ], pch = 16, asp = 1
  )
  lines(SlingshotDataSet(sce), lwd = 2, type = line_type)
  title(title)
  legend("bottomleft",
    legend = str_flatten_comma(included_clusters), title = "Clusters included:"
  )
  legend("topleft", legend = end_cluster, title = "End cluster:")
}

do_sling <- function(sce_sobj = sce_sobj,
                     end_cluster = NULL,
                     start_cluster = start_cluster) {
  sce <-
    slingshot(
      data = sce_sobj,
      clusterLabels = colData(sce_sobj)[["sailerx_clusters"]],
      # clusterLabels = sce_sobj %>% colData() %>% as_tibble() %>%
      # pull("sailerx_clusters"),
      reducedDim = "UMAP_SAILERX",
      start.clus = start_cluster,
      end.clus = end_cluster
    )
  return(sce)
}

do_heatmap <- function(features = VariableFeatures(sobj),
                       slingshot_df,
                       slot_to_use = "data",
                       show_row_names = TRUE,
                       cluster_features = FALSE,
                       supplied_data = NULL) {
  if (!is.null(supplied_data)) {
    mat_to_use <- supplied_data
  } else {
    mat_to_use <- slot(sobj@assays$RNA, slot_to_use)
  }
  subset_matrix <- mat_to_use[
    features,
    slingshot_df$barcodes
  ] %>%
    as.matrix()
  col_fun <- colorRamp2(
    c(0, quantile(subset_matrix, probs = c(0.99))[[1]]),
    c("darkblue", "green")
  )
  present_clusters <- unique(slingshot_df$sailerx_clusters)
  n_clusters <- (sobj[["sailerx_clusters"]] %>% unique() %>% dim())[1]
  sailerx_colours <- palette.colors(
    n = length(0:8), palette = "Classic Tableau"
  )
  names(sailerx_colours) <- as.character(0:8)
  cluster_annotation <- HeatmapAnnotation(
    df = slingshot_df["sailerx_clusters"],
    name = "sailerx_clusters",
    col = list(sailerx_clusters = sailerx_colours[present_clusters]),
    annotation_legend_param = list(
      legend_direction = "horizontal",
      title_position = "lefttop-rot"
    )
  )
  return(
    Heatmap(
      mat = subset_matrix,
      show_column_names = F,
      show_row_names = show_row_names,
      cluster_rows = cluster_features,
      cluster_columns = F,
      col = col_fun,
      column_title = "Cells ordered by pseudotime -->",
      row_title = "Features",
      heatmap_legend_param = list(
        title = "Expression",
        title_position = "lefttop-rot"
      ),
      border_gp = gpar(col = "black", lty = "solid"),
      column_order = slingshot_df$barcodes,
      top_annotation = cluster_annotation
    )
  )
}

correlation_plot <- function(end_cluster = NULL, quant = 0.90) {
  assert_that(!is.na(end_cluster))
  end_cluster <- as.character(end_cluster)
  quantil <- quantile(
    abs(correlations[[end_cluster]]),
    probs = quant, na.rm = T
  )
  ggplot(
    data.frame(
      "cor" = correlations[[end_cluster]],
      "high" = abs(correlations[[end_cluster]]) > quantil
    ),
    aes(seq_along(cor), y = cor, color = as.factor(high))
  ) +
    geom_line(
      mapping = aes(y = quantil),
      color = "gray20", linetype = "dashed"
    ) +
    geom_line(
      mapping = aes(y = -quantil),
      color = "gray20", linetype = "dashed"
    ) +
    geom_point() +
    scale_color_manual(
      values = c("darkgray", "maroon"),
      name = str_c("Over ", quant * 100, "% abs()"),
      breaks = c(F, T)
    ) +
    ggtitle(str_c(
      "Pseudotime correlation with features for\ntrajectory with end_cluster ",
      end_cluster
    )) +
    xlab("Variable features / HVGs") +
    ylab("Correlation with pseudotime") +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.line = element_line()
    ) +
    coord_fixed(ratio = 2e3)
}

correlation_heatmap <- function(end_cluster = 4,
                                quant = .98,
                                slot_to_use = "data",
                                cluster_features = FALSE,
                                supplied_data = NULL) {
  end_cluster <- as.character(end_cluster)
  quantil <- quantile(
    abs(correlations[[end_cluster]]),
    probs = quant, na.rm = T
  )
  highest_features <-
    VariableFeatures(sobj)[which(abs(correlations[[end_cluster]]) > quantil)]
  slingshot_df <- colData(sce[[traj_dict[[end_cluster]][["lvl1"]]]])
  slingshot_df[["slingshot"]] <- NULL
  slingshot_df <- slingshot_df %>% as.data.frame()
  slingshot_df <- as.data.frame(list(
    "barcodes" = rownames(slingshot_df),
    "sailerx_clusters" = slingshot_df$sailerx_clusters,
    "pseudotime" = slingPseudotime(
      sce[[traj_dict[[end_cluster]][["lvl1"]]]]
    )[, str_c("Lineage", traj_dict[[end_cluster]][["lvl2"]])]
    # "pseudotime" = slingshot_df[[str_c("slingPseudotime_",
    # traj_dict[[end_cluster]][["lvl2"]])]]
    # "pseudotime" = pseudotimes_of_interest[[end_cluster]]
  )) %>%
    arrange(pseudotime) %>%
    drop_na()
  do_heatmap(highest_features,
    slingshot_df = slingshot_df,
    slot_to_use = slot_to_use,
    cluster_features = cluster_features,
    supplied_data = supplied_data
  ) +
    labs(caption = end_cluster)
}

trajectory_plot <- function(end_cluster) {
  sds <- sce[[traj_dict[[as.character(end_cluster)]][["lvl1"]]]]
  nc <- 1
  pt <- slingPseudotime(sds)
  nms <- colnames(pt)
  nr <- ceiling(length(nms) / nc)
  pal <- viridis(100, end = 0.95)
  # par(mfrow = c(nr, nc))
  z <- traj_dict[[as.character(end_cluster)]][["lvl2"]] %>% as.numeric()
  i <- nms[z]
  colors <- pal[cut(pt[, i], breaks = 100)]
  plot(reducedDims(sce_all)$UMAP_SAILERX,
    pch = 16, col = rgb(red = .6, green = .6, blue = .6, alpha = .7), asp = 1
  )
  # plot(reducedDims(sds)$UMAP_SAILERX, pch = 16, main = i, col = rgb(red = .6,
  # green = .6, blue = .6, alpha = .7), asp = 1)
  points(reducedDims(sds)$UMAP_SAILERX,
    col = colors, pch = 16, cex = 0.5, asp = 1
  )
  lines(SlingshotDataSet(sds), lwd = 2, col = "black") # , type = 'lineages')
}
