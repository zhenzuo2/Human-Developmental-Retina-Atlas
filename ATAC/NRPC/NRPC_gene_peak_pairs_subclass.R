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

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.csv")

meta <- meta[meta$leiden %in% c("2", "13", "0", "1", "18", "9", "10", "4",
    "17", "7", "14"), ]

meta$subclass = "Unknown"
clusters = c("2", "13")
meta[meta$leiden %in% clusters, "subclass"] = "AC"
clusters = c("0", "1", "18")
meta[meta$leiden %in% clusters, "subclass"] = "HC"
clusters = c("9", "10")
meta[meta$leiden %in% clusters, "subclass"] = "BC"
clusters = c("4")
meta[meta$leiden %in% clusters, "subclass"] = "Rod"
clusters = c("17")
meta[meta$leiden %in% clusters, "subclass"] = "Cone"
clusters = c("7", "14")
meta[meta$leiden %in% clusters, "subclass"] = "RGC"

cells <- getCellNames(proj2)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X
common_cells <- intersect(cells, meta$X)

proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)
proj2$subclass <- meta[proj2$X2, "subclass"]

seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj2 <- addGeneExpressionMatrix(input = proj2, seRNA = seRNA, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "subclass")

proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "subclass",
    pathToMacs2 = findMacs2())

proj2 <- addPeakMatrix(proj2)

proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)

proj2 <- addPeak2GeneLinks(ArchRProj = proj2, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

proj2$Weeks = plyr::mapvalues(proj2$Days, from = c(70, 79, 87, 91, 100,
    103, 116, 120, 136, 137, 141, 142, 162, 165), to = c("PCW10", "PCW10",
    "PCW13", "PCW13", "PCW16", "PCW16", "PCW16", "PCW16", "PCW20", "PCW20",
    "PCW20", "PCW20", "PCW23", "PCW23"))

colp = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
names(colp) = c("PCW10", "PCW13", "PCW16", "PCW20", "PCW23")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Weeks", k = 6,
    palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/NRPC_gene_peak_pairs_Weeks.svg",
    width = 15, height = 8)
p
dev.off()

colp = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f", "#bcbd22")
names(colp) = c("AC", "BC", "Cone", "HC", "RGC", "Rod")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "subclass", k = 6,
    palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/NRPC_gene_peak_pairs_subclass.svg",
    width = 15, height = 8)
p
dev.off()

colp = c("#e6194B", "#000075")
names(colp) = c("Macula", "Peripheral")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Region", k = 6,
    palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/NRPC_gene_peak_pairs_Region.svg",
    width = 15, height = 8)
p
dev.off()

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-NRPC", load = FALSE)
