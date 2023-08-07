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

#time <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/AC_NRPC_NRPC_only_monocle3_pseudotime.csv",sep = " ")
#time$cut <- cut(time$x,3)
#proj2$cut <- time[stringi::stri_replace_last_fixed(rownames(proj2@cellColData),"#","_"),'cut']
#trajectory <- c("(-0.00589,1.96]","(1.96,3.93]","(3.93,5.9]")
#proj2 <- addTrajectory(
#    ArchRProj = proj2, 
#    name = "NRPC", 
#    groupBy = "cut",
#    trajectory = trajectory, 
#    embedding = "UMAP", 
#    force = TRUE
#)

trajectory <- c("PCW10", "PCW23")
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

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Weeks", embedding = "UMAP",size = 2)
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp.svg",
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
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=TRUE,varCutOff=0.6)
rownames(p1) <- str_extract(rownames(p1), "(?<=:)[A-Za-z0-9]+")
df <- read.csv("/storage/singlecell/zz4/fetal_snakemake/temp/highly_varible_gene_AC.csv")
p1 <- p1[rownames(p1) %in% df$X0,]
colnames(p1) <- NULL
svg(file="/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-AC_NRPC/Plots/Plot-MotifMatrix-Heatmaps.svg",height = 16,width = 12)
ht <- Heatmap(p1, name = "p1", row_order = 1:nrow(p1), column_order = 1:ncol(p1))
draw(ht)
dev.off()

trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE,groupEvery = 5)
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),returnMatrix=TRUE,varCutOff = 0.9)
rownames(p2) <- sub("chr\\d+:([^:]+)", "\\1", rownames(p2))
p2<- p2[rownames(p2) %in% rownames(p1),]
p2 <- p2[rownames(p1)[rownames(p1) %in% rownames(p2)],]
colnames(p2) <- NULL
svg(file="/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-AC_NRPC/Plots/Plot-GeneExpressionMatrix-Heatmaps.svg",height = 16,width = 12)
col_fun = colorRamp2(c(-2, 2), c("blue","yellow"))
ht <- Heatmap(p2, name = "p1", row_order = 1:nrow(p2), column_order = 1:ncol(p2),col = col_fun)
draw(ht)
dev.off()

assay(trajMM)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),varCutOff=0.6,labelTop = 20,grepExclude = removed)
#c("z:ONECUT1_295","z:ONECUT2_293","z:ONECUT3_296","z:NEUROD4_820","z:NEUROG3_43","z:NEUROD6_821","z:NEUROD2_73","z:NEUROD1_63","z:ISL1_398","z:PAX6_604","z:EBF1_67","z:BHLHE22_86","z:MEIS2_471","z:LHX5_406","z:LHX1_469","z:PTF1A_70","z:BARHL1_456")
plotPDF(p1, name = "Plot-MotifMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####
trajGEM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
p1 <- plotTrajectoryHeatmap(trajGEM, pal = paletteContinuous(set = "blueYellow"))
plotPDF(p1, name = "Plot-GeneExpressionMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####

trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"))
plotPDF(p1, name = "Plot-GeneScoreMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)
####
trajGSM <- getTrajectory(ArchRProj = proj2, name = "NRPC", useMatrix = "PeakMatrix", log2Norm = TRUE)
p1 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "greenBlue"))
plotPDF(p1, name = "Plot-PeakMatrix-Heatmaps.pdf", ArchRProj = proj2,
    addDOC = FALSE, width = 6, height = 8)

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
rownames(trajMM2) <- gsub(".*:(.*?)_.*", "\\1", rownames(trajMM2))
rownames(trajGEM2) <- gsub(".*:", "", rownames(trajGEM2))
ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = -1, rowOrder = rowOrder,labelTop = 100,labelMarkers = rownames(trajMM2))
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = -1, rowOrder = rowOrder,labelTop = 100,labelMarkers = rownames(trajMM2))
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-AC_NRPC/Plots/h1_h2.svg",width=18,height=15)
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

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-AC_NRPC", load = FALSE)