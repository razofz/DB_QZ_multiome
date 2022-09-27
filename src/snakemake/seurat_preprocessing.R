suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(stringr))
set.seed(snakemake@config[["seed"]])

###############
#  Load data  #
###############

data <- Read10X_h5(snakemake@input[["h5_10X"]])

sobj <- CreateSeuratObject(
  counts = data$`Gene Expression`,
  assay = "RNA"
)
sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = F)
seqlevelsStyle(annotation) <- "UCSC"
sobj[["ATAC"]] <- CreateChromatinAssay(
  counts = data$Peaks,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = snakemake@input[["fragments_file"]],
  annotation = annotation,
  verbose = F
)
DefaultAssay(sobj) <- "ATAC"
genome(sobj) <- "mm10"
sobj <- NucleosomeSignal(sobj, verbose = F)
sobj <- TSSEnrichment(sobj, verbose = F)
sobj <- subset(
  x = sobj,
  subset = nCount_ATAC < 1e5 &
    nCount_ATAC > 1000 &
    nCount_RNA < 5e4 &
    nCount_RNA > 1000 &
    nucleosome_signal < 3 &
    TSS.enrichment > 3
)

########################
#  Save Seurat object  #
########################

saveRDS(sobj, snakemake@output[["seurat_object"]])
