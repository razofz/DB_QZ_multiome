invisible(sapply(c(
  "stringr",
  "parallel",
  # "docstring",
  "patchwork",
  "ggplot2",
  "tidyr",
  "dplyr",
  # "tictoc",
  "readr",
  "assertthat",
  "tibble",
  # "circlize",
  # "slingshot",
  # "ComplexHeatmap",
  # "RColorBrewer",
  # "viridis",
  "Seurat",
  "EnsDb.Mmusculus.v79",
  "Signac"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

project_path <- Sys.getenv("PROJECT_PATH")
output_dir <- str_c(project_path, "/data/processed/notebooks/Immature/")
plot_dir <- str_c(output_dir, "plots/")
external_dir <- str_c(project_path, "/data/external/")
rds_path <- str_c(output_dir, "seurat_object_post_atac_processing.rds")
raw_dir <- str_c(project_path, "/data/raw/aggregate/2021_083_Qinyu_agg/outs/")
h5_path <- str_c(raw_dir, "filtered_feature_bc_matrix.h5")
fragments_path <- str_c(raw_dir, "atac_fragments.tsv.gz")
snakemake_analysis_path <- str_glue(
    "{project_path}/data/processed/snakemake/"
)

seed <- 12345
set.seed(seed)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# sailerx_colours <- palette.colors(n = 9, palette = "Paired")
# names(sailerx_colours) <- as.character(0:8)
