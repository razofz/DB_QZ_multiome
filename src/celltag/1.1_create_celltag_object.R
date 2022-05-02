lapply(list(
  "stringr",
  "CellTagR"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
})

set.seed(snakemake@config[["seed"]])
ct_version <- str_c("v", snakemake@config[["celltag_version"]])

local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cloud-project.org"
  options(repos = r)
})

bam_obj <- CellTagObject(
  object.name = str_c(snakemake@wildcards[["sample"]], "_bam"),
  fastq.bam.directory = snakemake@input[["bam"]]
)

# Extract the CellTag information
bam_obj <- CellTagExtraction(bam_obj, celltag.version = ct_version)
# Check the bam file result
print(head(bam_obj@bam.parse.rslt[[ct_version]]))

# Generate the sparse count matrix
bam_obj <- CellTagMatrixCount(
  celltag.obj = bam_obj,
  barcodes.file = snakemake@input[["barcodes_tsv"]]
  # need to add suffix -1 to barcodes file
  # okay, that is being done now, maybe I did smth manual there before
)
# Check the dimension of the raw count matrix
print(dim(bam_obj@raw.count))

# starcode part:
# Generating the collapsing file
bam_obj <- CellTagDataForCollapsing(
  celltag.obj = bam_obj, output.file = snakemake@output[["collapsing_file"]]
)

# save rds of bam_obj here, quit script, run starcode
saveRDS(bam_obj, snakemake@output[["bam_obj"]])
