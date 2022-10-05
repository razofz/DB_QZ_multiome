invisible(sapply(c(
  "TFBSTools",
  "stringr",
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

df <- read.table(snakemake@input[["annotations_bed"]], header = T)

################################################################################
#       Remove non-standard chromosomes, they break the adding of motifs       #
################################################################################

gr <- granges(sobj)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Mmusculus.UCSC.mm10)
seq_keep <- as.vector(seq_keep)
feat_keep <- GRangesToString(grange = gr[seq_keep])
sobj[["ATAC"]] <- subset(sobj[["ATAC"]], features = feat_keep)

# main_chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
# keep_peaks <- as.logical(seqnames(granges(sobj)) %in% main_chroms)
# sobj <- sobj[keep_peaks, ] # removes the RNA assay

gr <- granges(sobj)
seqlevels(gr) <- seqlevelsInUse(granges(sobj))


new_anno <- GRanges(
  ranges = IRanges(start = df$starts + 1, end = df$ends), df$seqname,
  gene_name = df$gene_name
)
genome(new_anno) <- unique(genome(Annotation(sobj[["ATAC"]])))

assay <- CreateChromatinAssay(
  fragments = Fragments(sobj[["ATAC"]]),
  motifs = Motifs(sobj[["ATAC"]]),
  seqinfo = seqinfo(sobj[["ATAC"]]),
  annotation = Annotation(sobj[["ATAC"]]),
  ranges = gr,
  data = sobj[["ATAC"]]@data,
  genome = "mm10"
)

sobj[["proximal"]] <- assay

# sobj@assays$ATAC@meta.features[peak_idx_proximal, "peak_type"] <- "proximal"

################################################################################
#              Split peaks into gene bodies, proximal and distal               #
################################################################################

# conflate all annotation regions for a gene, i.e. take the earliest start site
# and the latest end site and define that as the gene body

# classify peaks as "in gene body":
# should the whole peak be in the gene body?
# if not, how much of the peak should be in the gene body for the classification
# to happen? Should I consider if it "leaks" in the start or end? If in the
# start, might classify it as proximal/promoter instead.

promoters <- promoters(new_anno)
overlaps <- findOverlaps(gr, promoters)
peak_idx_proximal <- unique(overlaps@from)

# sobj@assays$ATAC@meta.features[peak_idx_proximal, "peak_type"] <- "proximal"
# sobj@assays$ATAC@meta.features[-peak_idx_proximal, "peak_type"] <-
#   "unclassified"


tmp <- rownames(sobj@assays$ATAC@meta.features[-peak_idx_proximal, ])
foo <- as.data.frame(t(data.frame(str_split(tmp, "-"))))
rownames(foo) <- 1:dim(foo)[1]
colnames(foo) <- c("seqnames", "start", "end")
unclassified <- GRanges(foo)
unclassified
in_gene_bodies <- findOverlaps(unclassified, new_anno)
peak_idx_gene_bodies <- unique(in_gene_bodies@from)
peak_idx_gene_bodies <- peak_idx_gene_bodies[-intersect(
  peak_idx_gene_bodies,
  peak_idx_proximal
)]

sobj@assays$ATAC@meta.features$peak_type <- "unclassified"
sobj@assays$ATAC@meta.features[peak_idx_gene_bodies, "peak_type"] <- "gene_body"
sobj@assays$ATAC@meta.features[peak_idx_proximal, "peak_type"] <- "proximal"
sobj@assays$ATAC@meta.features[
  sobj@assays$ATAC@meta.features$peak_type == "unclassified",
  "peak_type"
] <- "distal"



for (i in 45:50) {
  if (new_anno[overlaps@to[i]]$gene_name %in% rownames(sobj@assays$RNA)) {
    p <- CoveragePlot(sobj,
      region = new_anno[overlaps@to[i]]$gene_name,
      extend.upstream = 500, features = new_anno[overlaps@to[i]]$gene_name,
      extend.downstream = 500, peaks.group.by = "peak_type"
    )
  } else {
    p <- CoveragePlot(sobj,
      region = new_anno[overlaps@to[i]]$gene_name,
      extend.upstream = 500,
      extend.downstream = 500, peaks.group.by = "peak_type"
    )
  }
  p <- p +
    plot_annotation(title = str_c(
      "Visualising peaks around the ",
      new_anno[overlaps@to[i]]$gene_name,
      " body"
    ))
  ggsave(p,
    filename = str_c(
      "data/adhoc/ATAC/signac-motifs/",
      i, ".svg"
    ), device = "svg"
  )
}

feats <- VariableFeatures(sobj, assay = "RNA")
for (i in 10:20) {
  feat <- feats[i]
  p <- CoveragePlot(sobj,
    region = feat, features = feat, extend.upstream = 500,
    extend.downstream = 500, peaks.group.by = "peak_type",
    group.by = "origin"
  ) + plot_annotation(title = str_c(
    "Visualising peaks around the ",
    feat, " body"
  ))
  ggsave(p,
    filename = str_c(
      "data/adhoc/ATAC/signac-motifs/",
      "RNA_", i, ".svg"
    ), device = "svg"
  )
}



# classify peaks as "proximal":
# Region here is -2k:200 bp around TSS.
# But, same Q here: how much of the peak should be in this region?

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
  pfm = pfm
)

# system(str_c("touch ", snakemake@output[["seurat_object"]]))
saveRDS(snakemake@output[["seurat_object"]])
