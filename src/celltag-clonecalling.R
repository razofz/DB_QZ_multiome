suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(data.table))

#biddy_path <- snakemake@params[['biddy_path']]
biddy_path <- "src/BiddyetalWorkflow"
source(paste0(biddy_path, "/scripts/CellTagCloneCalling_Function.R"))

d13.mat <- as.data.frame(readRDS(snakemake@input[['matrix_rds']]))
# ./hf1.d15.v3.celltag.matrix.Rds
rownames(d13.mat) <- d13.mat$Cell.BC
d13.mat <- d13.mat[,-1]
d13.bin <- SingleCellDataBinarization(celltag.dat = d13.mat, 2)

d13.filt <- SingleCellDataWhitelist(celltag.dat = d13.bin, whitels.cell.tag.file = paste0(biddy_path, "/whitelist/V", snakemake@config[['celltag_version']], ".CellTag.Whitelist.csv"))
d13.filt <- MetricBasedFiltering(whitelisted.celltag.data = d13.filt, cutoff = 20, comparison = "less")
d13.filt <- MetricBasedFiltering(whitelisted.celltag.data = d13.filt, cutoff = 2, comparison = "greater")

d13.sim <- JaccardAnalysis(whitelisted.celltag.data = d13.filt, plot.corr = F, id = "d13", dir = dirname(snakemake@output[['clones_csv']]))
d13.clones <- CloneCalling(Jaccard.Matrix = d13.sim, output.dir = "", output.filename = snakemake@output[['clones_csv']], correlation.cutoff = 0.7)
#hf1.d15.v3.clones.csv

if (length(d13.clones) > 0) {
    write.csv(d13.clones[[2]], file=snakemake@output[['clones_size_csv']], row.names=F, quote=F)
} else {
    print(paste0('> No clones found.. Saving dummy outputs ',
		 snakemake@output[['clones_size_csv']], ' and ',
		 snakemake@output[['clones_csv']]))
    system(paste0('touch ', snakemake@output[['clones_size_csv']]))
    system(paste0('touch ', snakemake@output[['clones_csv']]))
}

