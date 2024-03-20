library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)
output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)

setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj <- loadArchRProject("Save-proj3")

p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE,
)
cmap  <- c('#9467bd', '#17becf', '#bcbd22', '#d62728', '#ff7f0e', '#e377c2', '#2ca02c', '#1f77b4', '#8c564b')
names(cmap) <- c('MG', 'Rod', 'BC', 'RGC', 'NRPC', 'Cone', 'HC', 'PRPC', 'AC')
p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "majorclass",palGroup = cmap)
svg("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure8/p2g_heatmap.svg")
p
dev.off()