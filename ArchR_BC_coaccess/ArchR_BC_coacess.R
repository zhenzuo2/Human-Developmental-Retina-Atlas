
library(ArchR)
library(Seurat)
input_path = "/storage/singlecell/zz4/fetal_bash/results/ArchR/"
input_rna_file = "/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
setwd(input_path)
projretina3 <- loadArchRProject("Save-projretina3")
cells <- rownames(projretina3@cellColData)
meta <- read.csv("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv")
rownames(meta) <- meta$X
cells = intersect(cells, stringi::stri_replace_last(meta$X, fixed = "_",
    "#"))
projretina4 <- subsetArchRProject(projretina3, cells = cells, force = TRUE)
seRNA <- readRDS("/storage/singlecell/zz4/fetal_bash/results/seRNA/seRNA.rds")
projretina4 <- addGeneExpressionMatrix(input = projretina4, seRNA = seRNA,
    force = TRUE)

# Co-accessibility with ArchR
projretina4 <- addCoAccessibility(ArchRProj = projretina4, reducedDims = "IterativeLSI")
cA <- getCoAccessibility(ArchRProj = projretina4, corCutOff = 0.5, resolution = 1,
    returnLoops = FALSE)
markerGenes <- c("OTX2", "VSX2", "ISL1", "PRDM1")
projretina4@cellColData$subclass <- meta[projretina4@cellColData$X, "subclass"]
projretina4@cellColData$majorclass <- meta[projretina4@cellColData$X, "majorclass"]
projretina4@cellColData$Days <- factor(meta[projretina4@cellColData$X,
    "Days"])
p <- plotBrowserTrack(ArchRProj = projretina4, groupBy = "subclass", geneSymbol = markerGenes,
    upstream = 50000, downstream = 50000, loops = getCoAccessibility(projretina4))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
    ArchRProj = projretina4, addDOC = FALSE, width = 5, height = 5)
p <- plotBrowserTrack(ArchRProj = projretina4, groupBy = "Days", geneSymbol = markerGenes,
    upstream = 50000, downstream = 50000, loops = getCoAccessibility(projretina4))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
    ArchRProj = projretina4, addDOC = FALSE, width = 5, height = 5)

# Peak2GeneLinkage with ArchR
projretina4 <- addPeak2GeneLinks(ArchRProj = projretina4, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix", )
p2g <- getPeak2GeneLinks(ArchRProj = projretina4, corCutOff = 0.45, resolution = 1,
    returnLoops = FALSE)
p <- plotPeak2GeneHeatmap(ArchRProj = projretina4, groupBy = "majorclass")
svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/ArchRSubset/Plots/Peak2GeneLinkage.svg")
p
dev.off()