invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat",
  "Signac"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])

DefaultAssay(sobj) <- "RNA"

################################################################################
#                                   Plotting                                   #
################################################################################

p <- DimPlot(sobj, group.by = "is_HSC") +
  coord_fixed()
ggsave(snakemake@output[["plot_is_HSC"]], plot = p, device = "svg")

p <- DimPlot(sobj, group.by = "is_earlydiff") + coord_fixed()
ggsave(snakemake@output[["plot_is_earlydiff"]], plot = p, device = "svg")
