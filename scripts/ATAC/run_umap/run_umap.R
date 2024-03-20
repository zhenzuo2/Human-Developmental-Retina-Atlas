output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(output_results_path)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj1 <- loadArchRProject("Save-proj2")

proj1 <- addIterativeLSI(ArchRProj = proj1, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)
proj1$Lab <- "Rui"
proj1@cellColData[proj1$Sample %in% c("sn_multiome_d59", "sn_multiome_d76p",
    "sn_multiome_d76c"), "Lab"] = "Tom"
proj1 <- addHarmony(ArchRProj = proj1, reducedDims = "IterativeLSI", name = "Harmony",
    groupBy = "Lab", corCutOff = 0.3, force = TRUE)

proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI", method = "Seurat",
    name = "Clusters", resolution = 0.8, force = TRUE)

proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "Harmony", name = "UMAPHarmony",
    nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)

p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "majorclass",
    embedding = "UMAPHarmony")
plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)

p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Lab",
    embedding = "UMAPHarmony")
plotPDF(p2, name = "Plot-UMAP-Sample-Clusters2.pdf", ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)

pathToMacs2 <- findMacs2()

proj1 <- addGroupCoverages(ArchRProj = proj1, groupBy = "majorclass")

proj1 <- addReproduciblePeakSet(ArchRProj = proj1, groupBy = "majorclass",
    pathToMacs2 = pathToMacs2)

proj1 <- addPeakMatrix(proj1)

markersPeaks <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "PeakMatrix",
    groupBy = "majorclass", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

proj1 <- addCoAccessibility(ArchRProj = proj1, reducedDims = "IterativeLSI")

cA <- getCoAccessibility(ArchRProj = proj1, corCutOff = 0.5, resolution = 1,
    returnLoops = FALSE)

proj1 <- addPeak2GeneLinks(ArchRProj = proj1, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

p2g <- getPeak2GeneLinks(ArchRProj = proj1, corCutOff = 0.45, resolution = 1,
    returnLoops = FALSE)

write.table(p2g, "/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g.csv",
    sep = ",", row.names = F, quote = F)
write.table(p2g@metadata$peakSet, "/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g_metadata_peakSet.csv",
    sep = ",", row.names = T, quote = F)
write.table(p2g@metadata$geneSet, "/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g_metadata_geneSet.csv",
    sep = ",", row.names = T, quote = F)

write.table(proj1@peakSet, "/storage/chentemp/zz4/adult_dev_compare/results/peaks/peakSet.csv",
    sep = ",", row.names = F, quote = F)

write.table(metadata(cA)$peakSet, "/storage/chentemp/zz4/adult_dev_compare/results/peaks/metadata_cA_peakSet.csv",
    sep = ",", row.names = F, quote = F)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj3", load = FALSE)