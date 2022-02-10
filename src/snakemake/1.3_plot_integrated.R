suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
set.seed(1234)

###############
#  Load data  #
###############

sobj <- readRDS(snakemake@input[["seurat_object"]])

##############
#  Plot RNA  #
##############

DefaultAssay(sobj) <- "RNA"

DimPlot(sobj, 



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
sobj <- RunHarmony(sobj, group.by.vars = "origin")

sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:30)
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:30) %>%
  FindClusters()

###############
#  ATAC part  #
###############

DefaultAssay(sobj) <- "ATAC"

sobj <- FindTopFeatures(sobj, min.cutoff = 5)
sobj <- RunTFIDF(sobj)
sobj <- RunSVD(sobj)
sobj <- RunHarmony(sobj, group.by.vars = "origin", reduction = "lsi", assay.use = "ATAC", project.dim = F)
sobj <- RunUMAP(sobj, dims = 2:30, reduction = "harmony")

########################
#  Save Seurat object  #
########################

saveRDS(sobj, snakemake@output[["seurat_object"]])
