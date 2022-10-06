invisible(sapply(c(
  "TFBSTools",
  "stringr",
  # "patchwork",
  # "ggplot2",
  "Seurat",
  "GenomicRanges",
  "BSgenome.Mmusculus.UCSC.mm10",
  "JASPAR2020",
  "Signac"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])
DefaultAssay(sobj) <- "ATAC"
Idents(sobj) <- "origin"

df <- read.delim(snakemake@input[["TSS_reference"]],
  sep = "\t", check.names = F
)

tss <- GRanges(
  seqnames = df$CHROM, IRanges(
    seqnames = df$CHROM, start = df$POS,
    end = df$POS, strand = df$STRAND
  ),
  symbol = df$Symbol, strand = df$STRAND, ID = df[["#Gene_IDs"]]
)
genome(tss) <- "mm10"

################################################################################
#       Remove non-standard chromosomes, they break the adding of motifs       #
################################################################################

gr <- granges(sobj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Mmusculus.UCSC.mm10)
seq_keep <- as.vector(seq_keep)
feat_keep <- GRangesToString(grange = gr[seq_keep])
sobj[["ATAC"]] <- subset(sobj[["ATAC"]], features = feat_keep)
gr <- gr[seq_keep]
genome(gr) <- "mm10"
seqlevels(gr) <- seqlevelsInUse(granges(sobj))
# sobj@assays$ATAC@meta.feature

################################################################################
#           Split peaks into proximal and distal (i.e. non-proximal)           #
################################################################################

promoters <- promoters(tss)
overlaps <- findOverlaps(gr, promoters)
peak_idx_proximal <- unique(overlaps@from)

sobj@assays$ATAC@meta.features$peak_type <- "unclassified"
sobj@assays$ATAC@meta.features[peak_idx_proximal, "peak_type"] <- "proximal"
sobj@assays$ATAC@meta.features[
  sobj@assays$ATAC@meta.features$peak_type == "unclassified",
  "peak_type"
] <- "distal"

proximal_peaks <-
  rownames(sobj@assays$ATAC@meta.features[
    sobj@assays$ATAC@meta.features$peak_type
    == "proximal",
  ])
distal_peaks <-
  rownames(sobj@assays$ATAC@meta.features[
    sobj@assays$ATAC@meta.features$peak_type
    == "distal",
  ])

sobj[["proximal"]] <- subset(sobj[["ATAC"]], features = proximal_peaks)
sobj[["distal"]] <- subset(sobj[["ATAC"]], features = distal_peaks)

#####################################################
#  Plot stuff (remove when running with snakemake)  #
#####################################################

# for (i in 9:15) {
#   should_save <- TRUE
#   if (is.null(LookupGeneCoords(
#     sobj,
#     tss[unique(overlaps@to)[i]]$symbol
#   ))) {
#     print(str_c(
#       "> ", tss[unique(overlaps@to)[i]]$symbol,
#       " was not found, could not generate plot."
#     ))
#     should_save <- FALSE
#   } else if (
#     tss[unique(overlaps@to)[i]]$symbol %in% rownames(sobj@assays$RNA)
#   ) {
#     p <- CoveragePlot(sobj,
#       region = tss[unique(overlaps@to)[i]]$symbol,
#       extend.upstream = 500, features = tss[unique(overlaps@to)[i]]$symbol,
#       extend.downstream = 500, peaks.group.by = "peak_type"
#     )
#   } else {
#     p <- CoveragePlot(sobj,
#       region = tss[unique(overlaps@to)[i]]$symbol,
#       extend.upstream = 500,
#       extend.downstream = 500, peaks.group.by = "peak_type"
#     )
#   }
#   if (should_save) {
#     p <- p +
#       plot_annotation(title = str_c(
#         "Visualising peaks around the ",
#         tss[unique(overlaps@to)[i]]$symbol,
#         " body"
#       ))
#     ggsave(p,
#       filename = str_c(
#         "data/adhoc/ATAC/signac-motifs/",
#         i, ".svg"
#       ), device = "svg"
#     )
#   }
# }

# feats <- VariableFeatures(sobj, assay = "RNA")
# for (i in 1:10) {
#   feat <- feats[i]
#   p <- CoveragePlot(sobj,
#     region = feat, features = feat, extend.upstream = 500,
#     extend.downstream = 500, peaks.group.by = "peak_type",
#     group.by = "origin"
#   ) + plot_annotation(title = str_c(
#     "Visualising peaks around the ",
#     feat, " body"
#   ))
#   ggsave(p,
#     filename = str_c(
#       "data/adhoc/ATAC/signac-motifs/",
#       "RNA_", i, ".svg"
#     ), device = "svg"
#   )
# }

################################################################################
#                                  Add motifs                                  #
################################################################################

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    all_versions = F
  )
)

sobj <- AddMotifs(
  object = sobj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "proximal",
  pfm = pfm
)

sobj <- AddMotifs(
  object = sobj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = "distal",
  pfm = pfm
)

################################################################################
#                           Compute motif activities                           #
################################################################################

sobj <- RunChromVAR(
    object = sobj,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    new.assay.name = str_c("proximal", "Chromvar"),
    assay = "proximal"
  )

sobj <- RunChromVAR(
    object = sobj,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    new.assay.name = str_c("distal", "Chromvar"),
    assay = "distal"
  )

saveRDS(sobj, snakemake@output[["seurat_object"]])
