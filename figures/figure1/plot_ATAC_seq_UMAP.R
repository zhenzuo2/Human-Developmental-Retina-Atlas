library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

proj2 <- loadArchRProject("Save-proj2")

proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

umap <- proj2@embeddings$UMAP$df
rownames(umap) <- stringi::stri_replace_last_fixed(rownames(umap),"#","_")
write.table(umap,"/storage/singlecell/zz4/fetal_snakemake/results/ArchR/atac_umap.csv")
saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-proj3", load = FALSE)