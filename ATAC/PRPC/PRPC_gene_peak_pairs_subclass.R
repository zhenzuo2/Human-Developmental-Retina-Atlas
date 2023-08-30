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

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_latent_time.csv")
meta$latent_time_cut <- cut(meta$latent_time,5)
cells <- getCellNames(proj2)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X
common_cells <- intersect(cells, meta$X)
proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

proj2$latent_time<- meta[common_cells,"latent_time"]
proj2$latent_time_cut<- meta[common_cells,"latent_time_cut"]

seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj2 <- addGeneExpressionMatrix(input = proj2, seRNA = seRNA, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "latent_time_cut")

proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "latent_time_cut",
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
names(colp) = c("(-0.001,0.2]","(0.2,0.4]", "(0.4,0.6]","(0.6,0.8]","(0.8,1]")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "latent_time_cut", k = 2,
    palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/PRPC_gene_peak_pairs_Weeks.svg",
    width = 10, height = 4)
p
dev.off()

p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "latent_time_cut", k = 2,
    palGroup = colp,returnMatrices=T)

colp = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f", "#bcbd22")
names(colp) = c("AC", "BC", "Cone", "HC", "RGC", "Rod")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "latent_time_cut", k = 2)
    #palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/NRPC_gene_peak_pairs_latent_time_cut.svg",
    width = 10, height = 8)
p
dev.off()

colp = c("#e6194B", "#000075")
names(colp) = c("Macula", "Peripheral")
p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "Region", k = 2,
    palGroup = colp)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/NRPC_gene_peak_pairs_Region.svg",
    width = 15, height = 8)
p
dev.off()


proj2 <- addCoAccessibility(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = proj2,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

markerGenes  <- c(
    "NFIA"
  )

p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "latent_time_cut", 
    geneSymbol = markerGenes, 
    upstream = 60000,
    downstream = 700000,
    loops = getPeak2GeneLinks(proj2)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 5, height = 5)

markerGenes  <- c("C1orf61")
p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "latent_time_cut", 
    geneSymbol = markerGenes, 
    loops = getPeak2GeneLinks(proj2,resolution = 1000)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 8, height = 5)

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-PRPC", load = FALSE)

proj2$latent_time_cut <- factor(proj2$latent_time_cut)
proj2 <- addCellColData(proj2,data = proj2$latent_time_cut,name = "latent_time_cut",cells = rownames(proj2@cellColData),force = T)

################################################################################################
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "PeakMatrix", 
    groupBy = "latent_time_cut",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

df <- markerList$'(0.8,1]'
df$peak <- paste(df$seqnames,df$start,df$end,sep="_")

p2g <- getPeak2GeneLinks(
    ArchRProj = proj2,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

p2g_links <- p2g

p2g_links$idxRNAName <- metadata(p2g)[[2]]$name[p2g_links$idxRNA]
temp <- as.data.frame(metadata(p2g)[[1]][p2g_links$idxATAC,])
p2g_links$idxATACName <- paste(temp$seqnames,temp$start,temp$end,sep="_")
p2g_links <- p2g_links[p2g_links$idxATACName %in% df$peak,]


#######
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

proj2 <- loadArchRProject("Save-PRPC")
cells <- rownames(proj2@cellColData[proj2$Sample %in% names(table(proj2$Sample)[table(proj2$Sample)>100]),])
proj2 <- loadArchRProject("Save-proj2")
proj3 <- subsetArchRProject(proj2,cells = cells, force = TRUE)

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_latent_time.csv")
meta$latent_time_cut <- cut(meta$latent_time,5)
cells <- getCellNames(proj3)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X

proj3$latent_time<- meta[cells,"latent_time"]
proj3$latent_time_cut<- meta[cells,"latent_time_cut"]


seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj3 <- addGeneExpressionMatrix(input = proj3, seRNA = seRNA, force = TRUE)

proj3 <- addGroupCoverages(ArchRProj = proj3, groupBy = "latent_time_cut")

proj3 <- addReproduciblePeakSet(ArchRProj = proj3, groupBy = "latent_time_cut",
    pathToMacs2 = findMacs2())

proj3 <- addPeakMatrix(proj3)

proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "cisbp", name = "Motif")
proj3 <- addBgdPeaks(proj3)
proj3 <- addDeviationsMatrix(
  ArchRProj = proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)

proj3 <- addCellColData(proj3,data = proj3$latent_time_cut,name = "latent_time_cut",cells = rownames(proj3@cellColData),force = T)
proj3 <- addIterativeLSI(ArchRProj = proj3, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)
proj3 <- addUMAP(
    ArchRProj = proj3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
trajectory = c("(0.2,0.4]","(0.8,1]")
proj3 <- addTrajectory(
    ArchRProj = proj3, 
    name = "PRPC", 
    groupBy = "latent_time_cut",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

trajMM  <- getTrajectory(ArchRProj = proj3, name = "PRPC", useMatrix = "MotifMatrix", log2Norm = FALSE,groupEvery = 5)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=TRUE)
rownames(p1) <- str_extract(rownames(p1), "(?<=:)[A-Za-z0-9]+")
write.csv(p1,"/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-PRPC/Plots/Plot-MotifMatrix-Heatmaps.csv")
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),labelTop=1000)
svg(file="/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-PRPC/Plots/Plot-MotifMatrix-Heatmaps.svg",height = 16,width = 12)
p1
dev.off()

saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-PRPC_motif", load = FALSE)


trajMM  <- getTrajectory(ArchRProj = proj3, name = "PRPC", useMatrix = "MotifMatrix", log2Norm = FALSE,groupEvery = 5)
trajGEM <- getTrajectory(ArchRProj = proj3, name = "PRPC", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE,groupEvery = 5)

corGSM_MM <- correlateTrajectories(trajGEM, trajMM,varCutOff1 = 0.7,varCutOff2 = 0.7,corCutOff = 0.4,log2Norm2 = FALSE)

trajGEM2 <- trajGEM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGEM2
assay(trajCombined, withDimnames=F) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))
rownames(trajMM2) <- gsub(".*:(.*?)_.*", "\\1", rownames(trajMM2))
rownames(trajGEM2) <- gsub(".*:", "", rownames(trajGEM2))
rownames(trajMM2) <- rownames(trajGEM2)
ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder,labelTop = 200)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder,force = TRUE,labelTop = 200)
svg("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/Save-PRPC_motif/Plots/h1_h2.svg",width=18,height=15)
ht1+ht2
dev.off()

ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder,labelTop = 200,returnMatrix=T)
svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/ht1.svg")
ComplexHeatmap::draw(ht1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()