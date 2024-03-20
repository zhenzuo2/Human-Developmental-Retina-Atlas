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

common_cells <- rownames(proj2@cellColData[proj2$majorclass=="Cone",])
proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj2 <- addGeneExpressionMatrix(input = proj2, seRNA = seRNA, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Region")

proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Region",
    pathToMacs2 = findMacs2())

proj2 <- addPeakMatrix(proj2)

proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)

proj2 <- addPeak2GeneLinks(ArchRProj = proj2, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "PeakMatrix", 
    groupBy = "Region",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "Region", 
    geneSymbol = c("CNGB3"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.5 & Log2FC >= 0.1", returnGR = TRUE)$Macula,
    upstream = 400000,
    downstream = 500000
)
plotPDF(p, name = "CNGB3", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "Region", 
    geneSymbol = c("PDE6H"),
    features =  getMarkers(markersPeaks, cutOff = "FDR < 0.5 & Log2FC >= 0.5", returnGR = TRUE)$Macula,
    upstream = 10000,
    downstream = 10000
)
plotPDF(p, name = "PDE6H", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Cone", load = FALSE)
