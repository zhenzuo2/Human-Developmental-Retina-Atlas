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
theme_ArchR(baseSize = 20)
# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-AC_NRPC")

trajectory <- c("PCW10", "PCW23")
proj2 <- addTrajectory(ArchRProj = proj2, name = "NRPC", groupBy = "Weeks",
    trajectory = trajectory, embedding = "UMAP", force = TRUE)

####
trajMM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "MotifMatrix",
    log2Norm = FALSE, groupEvery = 1)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, name = "Plot-MotifMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####
trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix",
    log2Norm = TRUE, groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),
    maxFeatures = 100, labelTop = 50)
plotPDF(p1, name = "Plot-GeneExpressionMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####
trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneScoreMatrix",
    log2Norm = TRUE, groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"),
    maxFeatures = 100, labelTop = 50)
plotPDF(p1, name = "Plot-GeneScoreMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####
trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "PeakMatrix",
    log2Norm = TRUE, groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "greenBlue"),
    maxFeatures = 100, labelTop = 50)
plotPDF(p1, name = "Plot-PeakMatrix-Heatmaps.pdf", ArchRProj = proj2, addDOC = FALSE,
    width = 6, height = 8)

proj3 <- subsetArchRProject(ArchRProj = proj2, cells = getCellNames(proj2)[proj2$Region ==
    "Peripheral"], outputDirectory = "ArchRSubset", force = TRUE)


trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix",
    log2Norm = TRUE, groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),
    maxFeatures = 1000, returnMatrix=T)
strings <-rownames(p1)
result <- sapply(strsplit(strings, ":"), function(x) x[2])
################################################################################################################################################################################################
markersPeaks <- getMarkerFeatures(ArchRProj = proj3, useMatrix = "PeakMatrix",
    groupBy = "Weeks", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
p <- plotBrowserTrack(ArchRProj = proj3, groupBy = "Weeks", geneSymbol = c("ONECUT1"),
    upstream = 20000, downstream = 20000, facetbaseSize = 20, baseSize = 12,
    scTileSize = 20)
plotPDF(p, name = "plotBrowserTrack_ONECUT1.pdf", ArchRProj = proj3, addDOC = FALSE,
    width = 6, height = 8)

p <- plotBrowserTrack(ArchRProj = proj3, groupBy = "Weeks", geneSymbol = c("PTF1A"),
    upstream = 5000, downstream = 5000, facetbaseSize = 20, baseSize = 12,
    scTileSize = 20)
plotPDF(p, name = "plotBrowserTrack_PTF1A.pdf", ArchRProj = proj3, addDOC = FALSE,
    width = 6, height = 8)

p <- plotBrowserTrack(ArchRProj = proj3, groupBy = "Weeks", geneSymbol = c("EBF3"),
    upstream = 50000, downstream = 50000, facetbaseSize = 20, baseSize = 12,
    scTileSize = 20)
plotPDF(p, name = "plotBrowserTrack_EBF3.pdf", ArchRProj = proj3, addDOC = FALSE,
    width = 6, height = 8)