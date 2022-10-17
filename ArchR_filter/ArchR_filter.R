args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
input_meta <- args[2]
output_results_path <- args[3]

library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(input_path)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

projretina1 <- loadArchRProject("Save-projretina1")

projretina2 <- filterDoublets(projretina1)

cells = getCellNames(projretina2)
meta.data <- read.table(input_meta)
meta.data$X2 <- stri_replace_last(rownames(meta.data), fixed = "_", "#")

common_cells <- intersect(cells, meta.data$X2)

projretina3 <- subsetArchRProject(projretina2, cells = common_cells, force = TRUE)
rownames(meta.data) <- meta.data$X2
for (col in colnames(meta.data)) {
    projretina3@cellColData[, col] <- meta.data[rownames(projretina3),
        col]
}
write.table(projretina3@cellColData, paste(output_results_path, "projretina3_cellColData.csv",
    sep = ""))
# add an Iterative LSI-based dimensionality reduction to an
# ArchRProject
projretina3 <- addIterativeLSI(ArchRProj = projretina3, useMatrix = "TileMatrix",
    name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(0.2),
        sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30,
    force = T)

# add Harmony Batch Corrected Reduced Dims to an ArchRProject
projretina3 <- addHarmony(ArchRProj = projretina3, reducedDims = "IterativeLSI",
    name = "Harmony", groupBy = "Sample", force = T)

# add cluster information to an ArchRProject
projretina3 <- addClusters(input = projretina3, reducedDims = "IterativeLSI",
    method = "Seurat", name = "Clusters", resolution = 0.8, force = T)

# add a UMAP embedding of a reduced dimensions object to an
# ArchRProject
projretina3 <- addUMAP(ArchRProj = projretina3, reducedDims = "IterativeLSI",
    name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)

addArchRThreads(threads = 1)
projretina3 <- addGroupCoverages(ArchRProj = projretina3, groupBy = "scpred_prediction",
    threads = 1)
addArchRThreads(threads = 20)

pathToMacs2 <- findMacs2()
projretina3 <- addReproduciblePeakSet(ArchRProj = projretina3, groupBy = "scpred_prediction",
    pathToMacs2 = pathToMacs2)
# Add a Peak Matrix to the ArrowFiles of an ArchRProject
projretina3 <- addPeakMatrix(projretina3)

saveArchRProject(ArchRProj = projretina2, outputDirectory = "Save-projretina2",
    load = FALSE)

saveArchRProject(ArchRProj = projretina3, outputDirectory = "Save-projretina3",
    load = FALSE)