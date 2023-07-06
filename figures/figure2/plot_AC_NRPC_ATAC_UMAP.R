output_results_path = "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
library(ggplot2)
library(viridis)
set.seed(0)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-AC_NRPC")

temp <- proj2@embeddings$UMAP$df
temp$Region <- proj2$Region
temp$Days <- proj2$Days
temp$Weeks <- as.factor(proj2$Weeks)
temp = temp[temp$Region=="Peripheral",]
myColors <- viridis(5)
names(myColors) <- levels(temp$Weeks)
p <- ggplot(temp, aes(x = temp$"IterativeLSI#UMAP_Dimension_1", y = temp$"IterativeLSI#UMAP_Dimension_2",
    colour = temp$Weeks,size =5)) + geom_point() + theme_bw() + scale_color_manual(values = myColors)+labs(y= "", x = "")  + labs(fill='') 
svg("/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_ATAC_UMAP_Peripheral.svg")
p
dev.off()

temp <- proj2@embeddings$UMAP$df
temp$Region <- proj2$Region
temp$Days <- proj2$Days
temp$Weeks <- as.factor(proj2$Weeks)
temp = temp[temp$Region=="Macula",]
myColors <- viridis(5)
names(myColors) <- levels(temp$Weeks)
p <- ggplot(temp, aes(x = temp$"IterativeLSI#UMAP_Dimension_1", y = temp$"IterativeLSI#UMAP_Dimension_2",
    colour = temp$Weeks,size =5)) + geom_point() + theme_bw() + scale_color_manual(values = myColors)+labs(y= "", x = "") + labs(fill='') 
svg("/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_ATAC_UMAP_Macula.svg")
p
dev.off()