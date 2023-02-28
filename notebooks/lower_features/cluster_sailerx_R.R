# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # Evaluate SAILERX results (in R)

# %%
source("helpers.R")

# %%
str_glue("- project_path: {project_path}")
str_glue("- output_dir: {output_dir}")
str_glue("- plot_dir: {plot_dir}")
str_glue("- seed: {seed}")

# %%
inputs <- c(
  "pre_SAILERX_rds" = str_glue("{output_dir}seurat_object_pre_SAILERX.rds"),
  "sailerx_embeddings" = str_glue("{output_dir}sailerx_embeddings.csv")
)
outputs <- c(
  "post_SAILERX_rds" = str_glue("{output_dir}seurat_object_post_SAILERX.rds")
)

# %%
"inputs:"
inputs
"outputs:"
outputs

# %% [markdown]
# ---

# %%
sobj <- readRDS(inputs["pre_SAILERX_rds"])

# %%
sobj

# %%
sobj@meta.data$pca_clusters <- sobj@meta.data$`RNA_snn_res.0.8`

# %%
sobj@meta.data$orig.ident <- "Diverse"

# %%
sobj[[]] %>% head

# %%
dir(inputs["pre_SAILERX_rds"] %>% dirname)

# %%
embeddings = read.csv(inputs["sailerx_embeddings"], row.names = "cells")#, col.names = as.character(1:51))

# %%
dim(embeddings)

# %%
colnames(embeddings) <- 1:50

# %%
head(embeddings)

# %% [markdown]
# ## Make new dimred for the seurat object

# %%
Reductions(sobj)

# %%
sobj[["SAILERX"]] <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(embeddings), key = "sailerx_", assay = "RNA", misc = list(notes = "Embeddings generated with SAILERX."))

# %%
Reductions(sobj)

# %%
DimPlot(sobj, reduction = "SAILERX", group.by = "pca_clusters") + coord_fixed()

# %%
DimPlot(sobj, reduction = "pca", group.by = "pca_clusters") + coord_fixed()

# %%
ProjectDim(sobj, reduction = "SAILERX")

# %%
ProjectDim(sobj, reduction = "pca")

# %%
slotNames(sobj@reductions$SAILERX)

# %%
dim(sobj@reductions$SAILERX@cell.embeddings)

# %% [markdown]
# ## RNA clustering + UMAP

# %%
DefaultAssay(sobj) <- "RNA"

# %%
sobj <- FindNeighbors(sobj, dims = 1:50, reduction = "SAILERX", graph.name = c("RNA_SAILERX_nn", "RNA_SAILERX_snn"))

# %%
Graphs(sobj)

# %%
sobj <- FindClusters(sobj, graph.name = "RNA_SAILERX_snn")

# %%
sobj[["sailerx_clusters"]] <- sobj[["seurat_clusters"]]

# %%
sobj <- RunUMAP(sobj, reduction = "SAILERX", dims = 1:50, reduction.name = "UMAP_SAILERX", reduction.key = "UMAPSAILERX_", return.model = T)

# %%
Reductions(sobj)

# %%
p = DimPlot(sobj, reduction = "UMAP_SAILERX", group.by = "sailerx_clusters", label = T, label.box = T, cols = "Paired") + coord_fixed()
p

# %%
DimPlot(sobj, reduction = "UMAP_SAILERX", group.by = "sailerx_clusters", label = T, label.box = T, cols = "Paired") + coord_fixed()

# %%
sobj@reductions %>% names

# %%
sobj@reductions$UMAP_SAILERX@cell.embeddings[,"UMAPSAILERX_2"] = -sobj@reductions$UMAP_SAILERX@cell.embeddings[,"UMAPSAILERX_2"]

# %%
sobj@reductions$UMAP_SAILERX@cell.embeddings %>% head

# %%
p = DimPlot(sobj, reduction = "UMAP_SAILERX", group.by = "sailerx_clusters", label = T, label.box = T, cols = "Paired") + coord_fixed()
p

# %%
n_clusters <- (sobj[["sailerx_clusters"]] %>% unique() %>% dim())[1]
sailerx_colours <- palette.colors(n = n_clusters, palette = "Classic Tableau")

# %%
p = DimPlot(sobj, reduction = "UMAP_SAILERX", group.by = "sailerx_clusters", label = T, label.box = T, cols = sailerx_colours) + coord_fixed()
p

# %%
ggsave(
  filename = str_c("umap_sailerx_RNA", ".svg"),
  path = plot_dir,
  plot = p,
  device = "svg",
  width = 7,
  height = 7,
  units = "in"
)

# %% [markdown]
# ## New UMAP w/ clusters separated

# %%
p <- DimPlot(sobj, reduction = "UMAP_SAILERX", split.by = "sailerx_clusters", group.by = "sailerx_clusters", label = T, label.box = T, cols = sailerx_colours, ncol = 3) + coord_fixed()
p

# %%
plots <- list()
paired_colours <- sailerx_colours
for (cl in 0:(n_clusters-1)) {
  plots[[as.character(cl)]] <- DimPlot(
    sobj,
    reduction = "UMAP_SAILERX",
    cells.highlight = WhichCells(
      sobj,
      expr = sailerx_clusters == cl
    ),
    cols.highlight = paired_colours[[cl+1]]
  ) +
    coord_fixed() +
    NoLegend() +
    # labs(
    #   title = cl,
    # ) +
    xlab("") +
    ylab("") +
    # NoAxes() #+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()#,
      # axis.ticks = element_blank()
    )
}

# %%
p <- (plots[["0"]] | plots[["1"]] | plots[["2"]]) /
  (plots[["3"]] | plots[["4"]] | plots[["5"]]) /
  (plots[["6"]] | plots[["7"]] | plots[["8"]])
p

# %%
ggsave(
  plot = p,
  filename = str_c(plot_dir, "umap_clusters_separately.svg"),
  device = "svg", height = 11, width = 11, units = "in"
)

# %% [markdown]
# ---
#
# ## Differentially expressed genes

# %%
p <- FeaturePlot(
  sobj,
  features = c(
    "nCount_RNA",
    "nCount_ATAC",
    "nFeature_RNA",
    "nFeature_ATAC"
  ),
  order = T,
  max.cutoff = "q99",
  reduction = "UMAP_SAILERX",
  coord.fixed = T
)
p

# %%
ggsave(
  filename = str_c("depth_stats", ".svg"),
  path = plot_dir,
  plot = p,
  device = "svg",
  width = 10,
  height = 10,
  units = "in"
)

# %%
markers_rna <- FindAllMarkers(sobj, assay = "RNA", features = VariableFeatures(sobj))

# %% tags=[]
markers_rna

# %% tags=[]
markers_rna

# %%
write.table(markers_rna, file = str_c(output_dir, "RNA_markers.tsv"), sep = "\t")

# %%
top_markers <- markers_rna %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top_markers

# %%
top_markers <- markers_rna %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
top_markers

# %%
only_markers = list()
for (i in unique(top_markers$cluster)) {
    only_markers[i] = c(top_markers[top_markers$cluster == i,"gene"])
}
only_markers

# %%
stopifnot(all(Idents(sobj) == sobj[[]][["sailerx_clusters"]] ))

# %%
plot_top_markers <- function(i) {
  n_clusters <- (sobj[["sailerx_clusters"]] %>% unique() %>% dim())[1]
  Idents(sobj) <- "sailerx_clusters"
  colours <- sailerx_colours
  if (typeof(i) != "character") {
    i <- as.character(i)
  }
  return(
    ((FeaturePlot(
      sobj,
      features = only_markers[[i]][1],
      reduction = "UMAP_SAILERX",
      max.cutoff = "q99"
    ) + coord_fixed()
      + xlab("")
      + ylab("")
    ) |
      (FeaturePlot(
        sobj,
        features = only_markers[[i]][2],
        reduction = "UMAP_SAILERX",
        max.cutoff = "q99"
      ) + coord_fixed()
        + xlab("")
        + ylab("")
      )) /
      (
        # (DimPlot(sobj, reduction = "UMAP_SAILERX", group.by = "sailerx_clusters", label = T, label.box = T, repel = T, cols = brewer.pal(n_clusters, "Paired")) + coord_fixed()) |
        (DimPlot(
          sobj,
          reduction = "UMAP_SAILERX",
          cells.highlight = WhichCells(sobj, expression = sailerx_clusters == i),
          cols.highlight = colours[(i %>% as.integer()) + 1]
        ) + coord_fixed() + NoLegend()) |
          (
            (VlnPlot(
              sobj,
              features = only_markers[[i]][1],
              cols = colours
            ) + NoLegend()
              + xlab("")
              + ylab("")
            ) /
              (VlnPlot(
                sobj,
                features = only_markers[[i]][2],
                cols = colours
              ) + NoLegend())
          )
      )
  ) + plot_annotation(title = str_c("Top markers for cluster ", i))
}

# %%
save_plot_top_markers <- function(i, modality = "RNA") {
  if (typeof(i) != "character") {
    i <- as.character(i)
  }
  ggsave(
    filename = str_c("top_markers_", modality, "_", i, ".svg"),
    path = plot_dir,
    plot = plot_top_markers(i),
    device = "svg",
    width = 10,
    height = 10,
    units = "in"
  )
}

# %%
for (i in unique(top_markers$cluster)) {
  save_plot_top_markers(i)
}

# %%
# save_plot_top_markers(0)

# %%
plot_top_markers(4)

# %%
plot_top_markers(2)

# %% tags=[]
plot_top_markers(8)

# %% [markdown]
# ## Save Seurat object as RDS

# %%
sobj

# %%
sobj %>% saveRDS(file = outputs[["post_SAILERX_rds"]], compress = F)

# %% [markdown]
# ## More FeaturePlots

# %%
FeaturePlot(sobj, features = "Meis1", max.cutoff = "q99", pt.size = .8) + coord_fixed()

# %%
FeaturePlot(sobj, features = "Plxdc2", max.cutoff = "q99", pt.size = .8) + coord_fixed()

# %%
FeaturePlot(sobj, features = c("Plxdc2", "nCount_ATAC"), ncol = 2, blend = T, max.cutoff = "q99", pt.size = .8, coord.fixed = T)

# %%
foo <- markers_rna %>%
dplyr::filter(cluster == 3) %>%
arrange(desc(avg_log2FC))
foo %>% head(n = 15)

# %%
?FeaturePlot

# %%
p = FeaturePlot(sobj, features = c("Mecom", "nFeature_RNA"), order = T, blend = T, max.cutoff = "q99", pt.size = .8, coord.fixed = T, combine = F)

p[1]
p[2]
p[3]
p[4]

# %%
p = FeaturePlot(
    sobj,
    features = c(foo$gene[1], "nFeature_RNA"),
    order = T,
    blend = T,
    max.cutoff = "q99",
    pt.size = 1.5,
    coord.fixed = T,
    # cols = c("gray", "chocolate4", "chartreuse"),
    combine = F
)

p[1]
p[3]

# %%
p = FeaturePlot(
    sobj,
    features = c(foo$gene[2], "nFeature_RNA"),
    order = T,
    blend = T,
    max.cutoff = "q99",
    pt.size = 1.5,
    coord.fixed = T,
    # cols = c("gray", "chocolate4", "chartreuse"),
    combine = F
)

p[1]
p[3]

# %%
p = FeaturePlot(sobj, features = c("Ngp", "nFeature_RNA"), order = T, blend = T, max.cutoff = "q99", pt.size = .8, coord.fixed = T, combine = F)
p[3]

# %%
p = FeaturePlot(sobj, features = c("Ngp", "nFeature_ATAC"), order = T, blend = T, max.cutoff = "q99", pt.size = .8, coord.fixed = T, combine = F)
p[3]

# %%
p = FeaturePlot(sobj, features = c("Ngp", "nCount_ATAC"), order = T, blend = T, max.cutoff = "q99", pt.size = .8, coord.fixed = T, combine = F)
p[3]

# %% [markdown]
# # Misc.

# %%
markers_rna %>% head

# %%
markers_rna %>%
    group_by(cluster) %>%
    summarise(
        `#` = n(),
        max_FC = max(avg_log2FC),
        min_FC = min(avg_log2FC),
        mean_FC = mean(avg_log2FC),
        median_FC = median(avg_log2FC),
        max_pct1 = max(pct.1),
        min_pct1 = min(pct.1),
    )

# %%
plot_hist <- function(i, feature = "avg_log2FC") {
  return(
    markers_rna %>%
      dplyr::filter(cluster == i) %>%
      ggplot(mapping = aes(x = .data[[feature]])) +
      geom_histogram(color = "white", bins = 40) +
      ggtitle(str_c("cluster ", i))
  )
}

# %%
p <- (plot_hist(0) | plot_hist(1) | plot_hist(2)) /
(plot_hist(3) | plot_hist(4) | plot_hist(5)) /
(plot_hist(6) | plot_hist(7) | plot_hist(8))
p

# %%
svg(filename = str_glue("{plot_dir}deg_histograms.svg"))
show(p)
dev.off()

# %%
plot_hist(2)

# %%
df <- sobj@reductions$UMAP_SAILERX@cell.embeddings %>%
  as.data.frame %>%
  rownames_to_column(var = "barcode") %>%
  as_tibble() %>%
  add_column("sailerx_clusters" = sobj[[]][["sailerx_clusters"]]) %>%
  add_column("is_cluster_i" = sobj[[]][["sailerx_clusters"]] == i)
df %>% head()

# %%
library(reshape2)

# %%
ggplot(
  data = df,
  mapping = aes(x = UMAPSAILERX_1, y = UMAPSAILERX_2, colour = is_cluster_i)
) +
  geom_point() +
  scale_colour_discrete(
  type = c("blue", "red")
      ) +
    # values = c("blue", "red")
  # scale_fill_manual(
  #   values = c("blue", "red"),
  #   limits = c(TRUE, FALSE),
  #   breaks = c(TRUE, FALSE)
  # ) +
  coord_fixed()
