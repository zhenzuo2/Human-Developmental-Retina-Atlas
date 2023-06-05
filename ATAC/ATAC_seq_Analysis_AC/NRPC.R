output_results_path = "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
set.seed(0)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj2")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC_AC.csv")

meta <- meta[meta$majorclass %in% c("NRPC"), ]

cells <- getCellNames(proj2)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X
common_cells <- intersect(cells, meta$X)

proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj2 <- addGeneExpressionMatrix(input = proj2, seRNA = seRNA, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Time")

proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Time",
    pathToMacs2 = findMacs2())

proj2 <- addPeakMatrix(proj2)

proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj2 <- addPeak2GeneLinks(ArchRProj = proj2, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

proj2$Weeks = plyr::mapvalues(proj2$Days, from = c(70, 79, 87, 91, 100,
    103, 116, 120, 136, 137, 141, 142, 162, 165), to = c("PCW10", "PCW10",
    "PCW13", "PCW13", "PCW16", "PCW16", "PCW16", "PCW16", "PCW20", "PCW20",
    "PCW20", "PCW20", "PCW23", "PCW23"))

trajectory <- c("PCW10", "PCW13", "PCW16", "PCW20", "PCW23")
proj2 <- addTrajectory(
    ArchRProj = proj2, 
    name = "NRPC", 
    groupBy = "Weeks",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)
p <- plotTrajectory(proj2, trajectory = "NRPC", colorBy = "cellColData", name = "NRPC")
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp.svg",
    width = 15, height = 8)
p[[1]]
dev.off()

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Weeks", embedding = "UMAP")
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp.svg",
    width = 15, height = 8)
p1
dev.off()



saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-AC_NRPC", load = FALSE)
