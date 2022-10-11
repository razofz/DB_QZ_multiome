lapply(list(
  "stringr",
  "CellTagR"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
})

set.seed(snakemake@config[["seed"]])
ct_version <- stringr::str_c("v", snakemake@config[["celltag_version"]])

bam_obj <- readRDS(snakemake@input[["bam_obj"]])

bam_obj <- JaccardAnalysis(bam_obj, fast = T)
# fast = T needed, otherwise error about proxy::simil not able to being typecast
# as dsTMatrix

# Call clones
bam_obj <- CloneCalling(celltag.obj = bam_obj, correlation.cutoff = 0.7)

# bam_obj@clone.composition[[ct_version]]
# bam_obj@clone.size.info[[ct_version]]

write.csv(
  bam_obj@clone.composition[[ct_version]],
  file = snakemake@output[["clones_csv"]]
)

# save it in the processed dir
saveRDS(bam_obj, snakemake@output[["bam_obj"]])

# seems like we can't use this, since we only have one time point read.
# Cf. their paper where they have three time points, with one version of CellTag
# used for each. That's hardcoded into the below function, so low probability on
# us using this function directly. Perhaps modified.

# bam_obj <- convertCellTagMatrix2LinkList(bam_obj)
