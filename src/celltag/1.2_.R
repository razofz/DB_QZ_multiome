lapply(list(
  "stringr",
  "CellTagR"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
})

set.seed(snakemake@config[["seed"]])
ct_version <- str_c("v", snakemake@config[["celltag_version"]])

bam_obj <- readRDS(snakemake@input[["bam_obj"]])

# Recount and generate collapsed matrix
bam_obj <- CellTagDataPostCollapsing(
  celltag.obj = bam_obj,
  collapsed.rslt.file = snakemake@input[["collapsing_result"]]
)
# Check the dimension of this collapsed count.
head(bam_obj@collapsed.count)

# Calling binarization
bam_obj <- SingleCellDataBinatization(bam_obj, 2)

svg(snakemake@output[["metric_plots_pre"]], width = 11, height = 11)
MetricPlots(bam_obj)
dev.off()

# whitelisting time, could do this specific to our data but then I would need
# the fastq file(s), so cannot do it now. Will have to use the (I guess) general
# whitelist included in the workflow repo.
# or, perhaps here I should use the new file Qinyu provided. Let's try both with
# snakemake

bam_obj <- SingleCellDataWhitelist(bam_obj, snakemake@input[["whitelist"]])

svg(snakemake@output[["metric_plots_post_whitelist"]], width = 11, height = 11)
MetricPlots(bam_obj)
dev.off()

bam_obj <- MetricBasedFiltering(bam_obj, 20, comparison = "less")
bam_obj <- MetricBasedFiltering(bam_obj, 1, comparison = "greater")

svg(
  snakemake@output[["metric_plots_post_metric_filtering"]],
  width = 11, height = 11
)
MetricPlots(bam_obj)
dev.off()

# name it post_filtering or smth
saveRDS(bam_obj, snakemake@output[["bam_obj"]])
