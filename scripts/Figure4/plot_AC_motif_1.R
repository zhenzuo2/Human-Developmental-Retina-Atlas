output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
set.seed(0)
setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj3")

meta <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/AC_w_NRPC.csv")

meta <- meta[meta$majorclass %in% c("NRPC"), ]

cells <- getCellNames(proj2)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X
common_cells <- intersect(cells, meta$X)

proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Time")

proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Time",
    pathToMacs2 = findMacs2())

proj2 <- addPeakMatrix(proj2)

proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30, force = TRUE)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

proj2 <- addPeak2GeneLinks(ArchRProj = proj2, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

proj2$Weeks = plyr::mapvalues(proj2$Days, from = c(59, 70, 76, 79, 87, 91, 100, 103, 116, 137, 141, 142, 162, 165), to = c("PCW8", "PCW10", "PCW10", "PCW10", "PCW13", "PCW13", "PCW15", "PCW15",
    "PCW15", "PCW19", "PCW19", "PCW19", "PCW23", "PCW23"))

trajectory <- c("PCW13", "PCW23")
proj2 <- addTrajectory(
    ArchRProj = proj2, 
    name = "NRPC", 
    groupBy = "Weeks",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)
p <- plotTrajectory(proj2, trajectory = "NRPC", colorBy = "cellColData", name = "NRPC",embedding="UMAP")
svg("/storage/chentemp/zz4/adult_dev_compare/temp/temp2.svg",
    width = 15, height = 8)
p[[1]]
dev.off()

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Weeks", embedding = "UMAP",size = 2)
svg("/storage/chentemp/zz4/adult_dev_compare/temp/temp.svg",
    width = 4, height = 4)
p1
dev.off()

proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
proj2 <- addBgdPeaks(proj2)
proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "Motif",
  force = TRUE
)
####
trajMM  <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "MotifMatrix", log2Norm = FALSE,groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=TRUE)
rownames(p1) <- str_extract(rownames(p1), "(?<=:)[A-Za-z0-9]+")
write.csv(p1,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-MotifMatrix-Heatmaps.csv")
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),labelTop=1000)
svg(file="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-MotifMatrix-Heatmaps.svg",height = 16,width = 12)
p1
dev.off()

trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE,groupEvery = 5)
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),returnMatrix=TRUE)
rownames(p2) <- sub("chr\\d+:([^:]+)", "\\1", rownames(p2))
write.csv(p2,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-GeneExpressionMatrix-Heatmaps.csv")
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),labelTop=50)
svg(file="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-GeneExpressionMatrix-Heatmaps.svg",height = 16,width = 12)
p2
dev.off()

trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE,groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, returnMatrix=TRUE)
write.csv(p1,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-GeneScoreMatrix-Heatmaps.csv")
p1 <- plotTrajectoryHeatmap(trajGSM,pal = paletteContinuous(set = "horizonExtra"),)
svg(file="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-GeneScoreMatrix-Heatmaps.svg",height = 16,width = 12)
p1
dev.off()
####
trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "PeakMatrix", log2Norm = TRUE,groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajGSM, returnMatrix=TRUE)
write.csv(p1,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-PeakMatrix-Heatmaps.csv")
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "greenBlue"))
svg(file="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/Plot-PeakMatrix-Heatmaps.svg",height = 16,width = 12)
p1
dev.off()

#############################################
trajMM  <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "MotifMatrix", log2Norm = FALSE,groupEvery = 5)
trajGEM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE,groupEvery = 5)

corGSM_MM <- correlateTrajectories(trajGEM, trajMM,varCutOff1 = 0.6,varCutOff2 = 0.6,corCutOff = 0.3,log2Norm2 = FALSE)

trajGEM2 <- trajGEM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGEM2
assay(trajCombined, withDimnames=F) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))

rownames(trajGEM2) <- gsub(".*:", "", rownames(trajGEM2))
mat <- trajGEM2@assays@data$smoothMat
rownames(mat) <- rownames(trajGEM2)
mat <- data.frame(mat[rowOrder,])
write.csv(mat,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/trajGEM2.csv")


rownames(trajMM2) <- gsub(".*:(.*?)_.*", "\\1", rownames(trajMM2))
mat <- trajMM2@assays@data$smoothMat
rownames(mat) <- rownames(trajMM2)
mat <- data.frame(mat[rowOrder,])
write.csv(mat,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/trajMM2.csv")

trajCombined <- trajGEM2
assay(trajCombined, withDimnames=F) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))
rownames(trajMM2) <- gsub(".*:(.*?)_.*", "\\1", rownames(trajMM2))
rownames(trajGEM2) <- gsub(".*:", "", rownames(trajGEM2))
ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = -1, rowOrder = rowOrder,labelTop = 100,labelMarkers = rownames(trajMM2))
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = -1, rowOrder = rowOrder,labelTop = 100,labelMarkers = rownames(trajMM2))
svg("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/Plots/h1_h2.svg",width=10,height=7)
ht1+ht2
dev.off()

###
motifPositions <- getPositions(proj2)
motifs <- c("ONECUT1","ONECUT2", "PTF1A")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Weeks")

seFoot <- getFootprints(
  ArchRProj = proj2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Weeks",
  flank = 200
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj2, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 10,
  flank = 200,
  pal = c("#4858A7", "#788FC8", "#F49B7C", "#B51F29")
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj2, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 10,
  flank = 200,
  pal = c("#4858A7", "#788FC8", "#F49B7C", "#B51F29")
)
###

markerGenes  <- c(
    "ONECUT2" 
  )
p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "Weeks", 
    geneSymbol = markerGenes, 
    upstream = 20000,
    downstream = 60000,
    loops = getCoAccessibility(proj2),
    sizes = c(12, 1.5, 3, 4),
    ylim = c(0, 1),
    pal = c("#4858A7", "#788FC8", "#F49B7C", "#B51F29")
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 5, height = 5)

proj2 <- addHarmony(ArchRProj = proj2, reducedDims = "IterativeLSI", name = "Harmony",
    groupBy = "Lab", corCutOff = 0.3, force = TRUE)
res <- proj2@embeddings$UMAPHarmony$df
res$majorclass <- proj2$majorclass
res$subclass <- proj2$subclass
res$Weeks <- proj2$Weeks
write.csv(res,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Save-AC_NRPC/ATAC_UMAPHarmony_cor.csv")


saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-AC_NRPC", load = FALSE)