suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(EnsDb.Mmusculus.v79))
suppressPackageStartupMessages(library(stringr))
set.seed(snakemake@config[["seed"]])


################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])

################################################################################
#                                      RNA part                                #
################################################################################

DefaultAssay(sobj) <- "RNA"

sobj <- NormalizeData(sobj, verbose = F)
sobj <- FindVariableFeatures(sobj, verbose = F)
sobj <- ScaleData(sobj, verbose = F)
sobj <- RunPCA(sobj, verbose = F)

###################################################
#  which dataset/sample did each cell come from?  #
###################################################

sobj[["is_Immature"]] <- str_detect(Cells(sobj), "-2")
origin <- sobj[["is_Immature"]]
origin <- unlist(origin)
origin[origin == T] <- "Immature"
origin[origin == F] <- "Diverse"
sobj[["origin"]] <- origin

############################
#  integrate with Harmony  #
############################

sobj <- RunHarmony(
  sobj,
  group.by.vars = "origin",
)
sobj@reductions$harmony_RNA <- sobj@reductions$harmony

sobj <- FindNeighbors(
  sobj,
  reduction = "harmony_RNA",
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
  reduction = "harmony_RNA",
  dims = 1:30,
  reduction.name = "UMAPharmonyRNA",
  reduction.key = "UMAPharmonyRNA_",
  nn.name = "NNharmonyRNA",
  assay = "RNA",
  umap.method = "umap-learn"
)

#############################################
#  Do the same for aggregated/unintegrated  #
#############################################

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
  reduction = "pca",
  dims = 1:30,
  reduction.name = "UMAPunintegratedRNA",
  reduction.key = "UMAPunintegratedRNA_",
  nn.name = "NNunintegratedRNA",
  assay = "RNA",
  umap.method = "umap-learn"
)

################################################################################
#                                  ATAC part                                   #
################################################################################

DefaultAssay(sobj) <- "ATAC"

sobj <- FindTopFeatures(sobj, min.cutoff = 5, verbose = F)
sobj <- RunTFIDF(sobj, verbose = F)
sobj <- RunSVD(sobj, verbose = F)
sobj <- RunHarmony(
  sobj,
  group.by.vars = "origin",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = F,
)
sobj@reductions$harmony_ATAC <- sobj@reductions$harmony

sobj <- FindNeighbors(
  sobj,
  reduction = "harmony_ATAC",
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
  reduction = "harmony_ATAC",
  dims = 2:30,
  reduction.name = "UMAPharmonyATAC",
  reduction.key = "UMAPharmonyATAC_",
  nn.name = "NNharmonyATAC",
  assay = "ATAC",
  umap.method = "umap-learn"
)

#############################################
#  Do the same for aggregated/unintegrated  #
#############################################

sobj <- FindNeighbors(
  sobj,
  reduction = "lsi",
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
  reduction = "lsi",
  dims = 2:30,
  reduction.name = "UMAPunintegratedATAC",
  reduction.key = "UMAPunintegratedATAC_",
  nn.name = "NNunintegratedATAC",
  assay = "ATAC",
  umap.method = "umap-learn"
)

################################################################################
#                              Save Seurat object                              #
################################################################################

saveRDS(sobj, snakemake@output[["seurat_object"]])
