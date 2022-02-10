suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(stringr))
set.seed(1234)

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
  fragments = snakemake@input[["fragments_file"]],
  annotation = annotation
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
##############
#  RNA part  #
##############

DefaultAssay(sobj) <- "RNA"

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj, verbose = FALSE)

# which dataset/sample did each cell come from?
sobj[["is_Immature"]] <- str_detect(Cells(sobj), "-2")
origin <- sobj[["is_Immature"]]
origin <- unlist(origin)
origin[origin == T] <- "Immature"
origin[origin == F] <- "Diverse"
sobj[["origin"]] <- origin

# integrate with Harmony
sobj <- RunHarmony(
  sobj,
  group.by.vars = "origin",
  reduction.save <- "harmonyRNA"
)

sobj <- FindNeighbors(
  sobj,
  reduction = "harmonyRNA",
  dims = 1:30,
  assay = "RNA",
  graph.name = c("NNharmonyRNA", "SNNharmonyRNA")
)
sobj <- FindClusters(
  sobj,
  graph.name = "SNNharmonyRNA"
)
sobj[["clusters_RNA_harmony"]] <- Idents(sobj)
sobj <- RunUMAP(
  sobj,
  reduction = "harmonyRNA",
  dims = 1:30,
  reduction.key = "UMAPharmonyRNA_",
  nn.name = "NNharmonyRNA",
  graph = "SNNharmonyRNA",
  assay = "RNA",
  umap.method = "umap-learn"
)

# Do the same for aggregated/unintegrated
sobj <- FindNeighbors(
  sobj,
  reduction = "pca",
  dims = 1:30,
  assay = "RNA",
  graph.name = c("NNunintegratedRNA", "SNNunintegratedRNA")
)
sobj <- FindClusters(
  sobj,
  graph.name = "SNNunintegratedRNA"
)
sobj[["clusters_RNA_unintegrated"]] <- Idents(sobj)
sobj <- RunUMAP(
  sobj,
  reduction = "unintegratedRNA",
  dims = 1:30,
  reduction.key = "UMAPunintegratedRNA_",
  nn.name = "NNunintegratedRNA",
  graph = "SNNunintegratedRNA",
  assay = "RNA",
  umap.method = "umap-learn"
)

###############
#  ATAC part  #
###############

DefaultAssay(sobj) <- "ATAC"

sobj <- FindTopFeatures(sobj, min.cutoff = 5)
sobj <- RunTFIDF(sobj)
sobj <- RunSVD(sobj)
sobj <- RunHarmony(
  sobj,
  group.by.vars = "origin",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = F,
  reduction.save = "harmonyATAC"
)

sobj <- FindNeighbors(
  sobj,
  reduction = "harmonyATAC",
  dims = 1:30,
  assay = "ATAC",
  graph.name = c("NNharmonyATAC", "SNNharmonyATAC")
)
sobj <- FindClusters(
  sobj,
  graph.name = "SNNharmonyATAC"
)
sobj[["clusters_ATAC_harmony"]] <- Idents(sobj)
sobj <- RunUMAP(
  sobj,
  reduction = "harmonyATAC",
  dims = 2:30,
  reduction.key = "UMAPharmonyATAC_",
  nn.name = "NNharmonyATAC",
  graph = "SNNharmonyATAC",
  assay = "ATAC",
  umap.method = "umap-learn"
)

# Do the same for aggregated/unintegrated
sobj <- FindNeighbors(
  sobj,
  reduction = "unintegratedATAC",
  dims = 1:30,
  assay = "ATAC",
  graph.name = c("NNunintegratedATAC", "SNNunintegratedATAC")
)
sobj <- FindClusters(
  sobj,
  graph.name = "SNNunintegratedATAC"
)
sobj[["clusters_ATAC_unintegrated"]] <- Idents(sobj)
sobj <- RunUMAP(
  sobj,
  reduction = "unintegratedATAC",
  dims = 2:30,
  reduction.key = "UMAPunintegratedATAC_",
  nn.name = "NNunintegratedATAC",
  graph = "SNNunintegratedATAC",
  assay = "ATAC",
  umap.method = "umap-learn"
)

########################
#  Save Seurat object  #
########################

saveRDS(sobj, snakemake@output[["seurat_object"]])
