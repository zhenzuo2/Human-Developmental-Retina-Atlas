output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(output_results_path)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 1)

proj1 <- loadArchRProject("Save-proj3")

markersPeaks <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "PeakMatrix",
    groupBy = "majorclass", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

proj1@cellColData$majorclass <- factor(proj1$majorclass, levels = c("AC",
    "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"))
p <- plotBrowserTrack(ArchRProj = proj1, groupBy = "majorclass", geneSymbol = c("PDE6H"),
    features = getMarkers(markersPeaks, cutOff = "FDR <= 1.0 & Log2FC >= 0.5",
        returnGR = TRUE)["Cone"], upstream = 10000, downstream = 20000,
    pal = c("#8c564b", "#bcbd22", "#e377c2", "#2ca02c", "#9467bd", "#ff7f0e",
        "#1f77b4", "#d62728", "#17becf"))
plotPDF(p, name = "Plot-Tracks-With-Features-PDE6H", width = 5, height = 5,
    ArchRProj = proj1, addDOC = FALSE)

##########################################################################################

cells <- rownames(proj1@cellColData)[(proj1@cellColData$majorclass == "Cone") &
    (proj1@cellColData$Region != "Whole Eye")]

proj2 <- subsetArchRProject(proj1, cells = cells, force = TRUE)

markersPeaks <- getMarkerFeatures(ArchRProj = proj2, useMatrix = "PeakMatrix",
    groupBy = "Region", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")

p <- plotBrowserTrack(ArchRProj = proj2, groupBy = "Region", geneSymbol = c("PDE6H"),
    features = getMarkers(markersPeaks, cutOff = "FDR <= 1.0 & Log2FC >= 0.5",
        returnGR = TRUE)["Macula"], upstream = 10000, downstream = 20000,
    pal = c("#1f77b4", "#ff7f0e"))
plotPDF(p, name = "Plot-Tracks-With-Features-PDE6H-Region", width = 5,
    height = 5, ArchRProj = proj1, addDOC = FALSE)

