library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)
output_results_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)

setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj <- loadArchRProject("Save-proj3")

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix",
    groupBy = "majorclass", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")

file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/markersPeaks.rds"
saveRDS(
  markersPeaks,
  file_path
)

heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    labelRows = FALSE, pal = gplots::bluered(100))

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6,
    ArchRProj = proj, addDOC = FALSE)