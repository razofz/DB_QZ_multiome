invisible(sapply(c(
  "stringr",
  "Signac",
  "GenomicRanges",
  "EnsDb.Mmusculus.v79",
  "future.apply"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = F)
seqlevelsStyle(annotation) <- "UCSC"

genes <- unique(annotation$gene_name)
genes <- genes[str_length(genes) > 2]

plan(multisession, workers = 16)

annotations <- unlist(GRangesList(future_lapply(genes, FUN = function(x) {
  obj <- annotation[annotation$gene_name == x]
  if (length(obj) > 0) {
    ranges <- GRanges(str_c(
      as.character(unique(seqnames(obj))),
      ":", min(start(obj)), "-",
      max(end(obj))
    ),
    gene_name = unique(obj$gene_name),
    # gene_id = unique(obj$gene_id)
    )
    # genome(ranges) <- unique(genome(obj))
    # seqlengths(ranges) <- seqlengths(obj)[as.character(seqnames(ranges))]
    return(ranges)
  }
})))
annotations

df <- data.frame(
  seqname = seqnames(annotations), starts = start(annotations) - 1,
  ends = end(annotations), names = c(rep(".", length(annotations))),
  scores = c(rep(".", length(annotations))),
  strands = strand(annotations), gene_name = annotations$gene_name
)

write.table(df,
  file = snakemake@output[["annotations_bed"]],
  quote = F, sep = "\t", row.names = F,
  col.names = T
)
