library(ArchR)
library(ggplot2)
library(ggforce)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

input_path = "/storage/singlecell/zz4/fetal_bash/results/ArchR/"
set.seed(1)
setwd(input_path)

markerGenes <- c("FOXN4", "PROX1", "WIF1", "PTN", "ABI3BP", "GPX3", "NPVF",
    "TF", "LIMD1", "NFIX")

PRPC <- loadArchRProject("PRPC")
PRPC$Days <- as.character(PRPC$Days)
PRPC$Time_meta <- mapvalues(PRPC$Days, from = c(70, 79, 87, 91, 100, 103,
    116, 136, 137, 141, 142, 162, 165), to = c("FW10", "FW10", "FW13",
    "FW13", "FW13", "FW16", "FW16", "FW19", "FW19", "FW19", "FW19", "FW23",
    "FW23"), warn_missing = TRUE)

# Add pseudotime
Time <- read.csv("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv")
Time$X <- stringi::stri_replace_last(Time$X, fixed = "_", "#")
rownames(Time) <- Time$X
common_cells <- intersect(PRPC$cellNames, Time$X)
Time <- Time[common_cells, ]
PRPC <- subsetArchRProject(PRPC, cells = common_cells, force = TRUE)

Time[Time$latent_time <= quantile(Time$latent_time, 0.25), "Time_meta"] = "1"
Time[(Time$latent_time >= quantile(Time$latent_time, 0.25)) & (Time$latent_time <=
    quantile(Time$latent_time, 0.5)), "Time_meta"] = "2"
Time[(Time$latent_time >= quantile(Time$latent_time, 0.5)) & (Time$latent_time <=
    quantile(Time$latent_time, 0.75)), "Time_meta"] = "3"
Time[Time$latent_time >= quantile(Time$latent_time, 0.75), "Time_meta"] = "4"
PRPC@cellColData$Time_meta <- Time[rownames(PRPC@cellColData), "Time_meta"]
PRPC@cellColData$latent_time <- Time[rownames(PRPC@cellColData), "latent_time"]

PRPC_early <- subsetArchRProject(PRPC, cells = rownames(Time)[Time$Time_meta %in%
    c("1", "2")], force = TRUE, outputDirectory = "PRPC_early")
PRPC_late <- subsetArchRProject(PRPC, cells = rownames(Time)[Time$Time_meta %in%
    c("3", "4")], force = TRUE, outputDirectory = "PRPC_late")

# Plot UMAP by time
UMAP <- PRPC@embeddings$UMAP$df
UMAP$color <- as.numeric(PRPC$latent_time)
colnames(UMAP) <- c("x", "y", "color")

p <- ggplot(UMAP, aes(x = -x, y = y, colour = color)) + geom_point() +
    scale_color_viridis_c() + theme_bw() + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.background = element_blank())

svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/PRPC/Plots/IterativeLSI#UMAP.svg")
p
dev.off()

# addCoAccessibility
PRPC <- addCoAccessibility(ArchRProj = PRPC, reducedDims = "IterativeLSI")

cA <- getCoAccessibility(ArchRProj = PRPC, corCutOff = 0.5, resolution = 1,
    returnLoops = FALSE)
markerGenes = "NFIB"
p <- plotBrowserTrack(ArchRProj = PRPC, groupBy = "Time_meta", geneSymbol = markerGenes,
    sizes = c(10, 1.5, 3, 4), upstream = 4e+05, downstream = 1e+05, loops = getCoAccessibility(PRPC),
    , ylim = c(0, 1))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
    ArchRProj = PRPC, addDOC = FALSE, width = 5, height = 5)

# addPeak2GeneLinks
PRPC_early <- addPeak2GeneLinks(ArchRProj = PRPC_early, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = PRPC_early, corCutOff = 0.45, resolution = 1,
    returnLoops = FALSE)
p2gearly <- paste(p2g$idxATAC, p2g$idxRNA, sep = "_")
p <- plotBrowserTrack(ArchRProj = PRPC_early, groupBy = "Time_meta", geneSymbol = markerGenes,
    upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(PRPC_early))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
    ArchRProj = PRPC_early, addDOC = FALSE, width = 5, height = 5)
p <- plotPeak2GeneHeatmap(ArchRProj = PRPC_early, groupBy = "Time_meta")
svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/PRPC_early/Plots/Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.svg")
p
dev.off()

PRPC_late <- addPeak2GeneLinks(ArchRProj = PRPC_late, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = PRPC_late, corCutOff = 0.45, resolution = 1,
    returnLoops = FALSE)
p2glate <- paste(p2g$idxATAC, p2g$idxRNA, sep = "_")
p <- plotBrowserTrack(ArchRProj = PRPC_late, groupBy = "Time_meta", geneSymbol = markerGenes,
    upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(PRPC_late))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf",
    ArchRProj = PRPC_late, addDOC = FALSE, width = 5, height = 5)
p <- plotPeak2GeneHeatmap(ArchRProj = PRPC_late, groupBy = "Time_meta")
svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/PRPC_late/Plots/Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.svg")
p
dev.off()

# Define the two vectors of strings
vec1 <- p2gearly
vec2 <- p2glate

set1 <- setNames(as.list(vec1), vec1)
set2 <- setNames(as.list(vec2), vec2)

# Create Venn diagram
venn <- venn.diagram(list(vec1 = set1, vec2 = set2), filename = "/storage/singlecell/zz4/fetal_bash/results/ArchR/PRPC_late/Plots/Venn_diagram.png",
    fill = c("cornflowerblue", "darkorange"), alpha = c(0.5, 0.5), cat.fontface = "bold",
    cat.fontsize = 16, cat.dist = c(0.03, -0.03), cat.default.pos = "outer",
    cat.pos = c(0, 0), scaled = TRUE,category.names = c("Early" , "Late"))

# markersPeaks
markersPeaks <- getMarkerFeatures(ArchRProj = PRPC, useMatrix = "PeakMatrix",
    groupBy = "Time_meta", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6,
    ArchRProj = PRPC, addDOC = FALSE)

# addTrajectory
trajectory <- c("1", "2", "3", "4")
PRPC <- addTrajectory(ArchRProj = PRPC, name = "MyeloidU", groupBy = "Time_meta",
    trajectory = trajectory, embedding = "UMAP", force = TRUE)
head(PRPC$MyeloidU[!is.na(PRPC$MyeloidU)])
p <- plotTrajectory(PRPC, trajectory = "MyeloidU", colorBy = "cellColData",
    name = "MyeloidU")
plotPDF(p, name = "Plot-MyeloidU-Traj-UMAP.pdf", ArchRProj = PRPC, addDOC = FALSE,
    width = 5, height = 5)
PRPC <- addPeakMatrix(PRPC, force = TRUE)
PRPC <- addBgdPeaks(PRPC)
PRPC <- addDeviationsMatrix(ArchRProj = PRPC, peakAnnotation = "Motif",
    force = TRUE)
trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "MotifMatrix",
    log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, name = "Plot-MyeloidU-Traj-Heatmaps.pdf", ArchRProj = projHeme5,
    addDOC = FALSE, width = 6, height = 8)

# getMarkerFeatures (GeneScoreMatrix)
markersGS <- getMarkerFeatures(ArchRProj = PRPC, useMatrix = "GeneScoreMatrix",
    groupBy = "Time_meta", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
labelMarkers = c("CREB5", "E2F2", "E2F7", "E2F8", "EGR1", "ESRRG", "FOS",
    "FOXM1", "FOXN4", "FOXP2", "GLIS3", "HES6", "HMGB2", "ID2", "MECOM",
    "MXD3", "MYBL1", "NFIA", "NFIB", "NFIX", "NPAS3", "NR4A1", "PROX1",
    "RARB", "SOX5", "SOX6", "TCF4", "TOX", "TOX3", "ZIC1", "ZNF367", "ZNF730",
    "ABI3BP", "ALDH1A2", "BLM", "BRCA1", "BRIP1", "CENPK", "CGNL1", "CLSPN",
    "CNTN5", "DNA2", "DST", "EZH2", "FSTL5", "GRIA1", "KDM5B", "LRRC4C",
    "MMS22L", "NCALD", "ROBO2", "TRPM3", "UHRF1", "VRK1", "WDHD1", "XRCC2",
    "ASF1B", "ATAD2", "CHAF1A", "DTL", "ESCO2", "FANCA", "FANCD2", "KIF15",
    "MND1", "NCAPD3", "RRM2", "ZGRF1", "ADGRB3", "APOLD1", "AURKB", "BIRC5",
    "BUB1", "BUB1B", "CDCA3", "CDK1", "CENPE", "CENPF", "CKAP2", "COL25A1",
    "DEPDC1B", "ECT2", "EFNA5", "FAM110B", "H2AFX", "HHIP", "HIST1H1B",
    "KIF11", "KIF18B", "KIFC1", "KNL1", "LMNB1", "MAPK10", "MKI67", "NCAPG",
    "PBK", "PTMA", "RACGAP1", "RGS16", "RGS3", "RHBDL3", "RTKN2", "SGO2",
    "SMOC1", "SRRM4", "TOP2A", "TPX2", "TTK", "TUBA1B", "CADM1", "GRID2",
    "KIRREL3", "LINC00511", "NRG1", "PLPPR1", "RALYL", "RELN", "SPON1",
    "WNT5B", "FIGN", "FLT1", "NAV2", "PCDH7", "WIF1", "AGBL4", "AHI1",
    "CABLES1", "DCDC1", "DNAH11", "EPHA3", "FAP", "GALNT13", "GPC6", "KIAA1217",
    "MAN1A1", "MYO16", "OPCML", "PALLD", "PDLIM5", "PTPRM", "RMST", "SRGAP1",
    "APBB2", "ARHGAP32", "ARHGEF3", "ASPH", "ATP2B2", "BTG2", "C1orf61",
    "CD44", "CDON", "CHL1", "CNTN1", "COBL", "COL23A1", "COL9A1", "CRYAB",
    "CYP1B1", "DCC", "DCLK1", "DGKB", "DIAPH3", "DNER", "EEPD1", "ELL2",
    "ERC2", "F3", "FAM160A1", "FRMD5", "GABRG3", "GFRA1", "GPX3", "GRB10",
    "GRIK2", "KAZN", "KCNA2", "LGR5", "LIMD1", "LINC01619", "LRP2", "LSAMP",
    "MAP2", "MAPK4", "MCAM", "MIR217HG", "MTUS2", "NAV3", "NCKAP5", "NPVF",
    "NRXN3", "OSBPL10", "PAG1", "PCSK2", "PLCE1", "PLEKHG1", "PRKCA", "PTN",
    "QDPR", "RASD1", "RASGRP1", "RGR", "RLBP1", "SAMD5", "SCARB1", "SDK1",
    "SEL1L3", "SEMA5B", "SLC15A3", "SLC1A3", "SLC35F1", "SLC6A11", "SPOCK1",
    "SPTLC3", "ST6GAL2", "STK39", "TENM3", "TF", "TNIK", "UTRN", "VIPR2",
    "ZFHX4", "ADGRV1", "FBXO5", "IQGAP3", "MELK", "NCAPG2", "NCAPH", "NDC80",
    "POC1A", "RBL1", "SKA3", "SPC24", "WDR62", "WDR76", "FANCI", "ADD3",
    "ARHGAP19", "ASPM", "CDC25C", "CDCA2", "CENPU", "CEP112", "CEP128",
    "CIT", "CKAP5", "GINS1", "KIF14", "KIF23", "KIF24", "LINC01572", "NUCKS1",
    "PARPBP", "PIF1", "PRR11", "RASGEF1B", "SASH1", "SGO1", "SMC4", "SPC25",
    "STIL", "TACC3", "TMPO", "TRIM59", "UBE2T", "C21orf58", "CKAP2L", "CMC2",
    "DLGAP5", "NCAPD2", "PRC1", "CPA6", "HKDC1", "MDK", "MEG8", "MYO5B",
    "P3H2", "SERINC5", "ST3GAL1", "ST5", "ST8SIA5", "STAC", "CDH23", "CDH6",
    "NHSL1", "ANLN", "CHEK1", "HELLS", "MSH5", "NUSAP1", "ALDH1A1", "CALM1",
    "CHSY3", "CPAMD8", "EGFR", "EYA2", "FRMD4B", "LAMA1", "MIS18BP1", "PDE4D",
    "PHLPP1", "PTPRO", "RBFOX1", "SHROOM3", "SLC6A1", "SLC8A1", "TRIM2",
    "BARD1", "CCDC18", "CCNF", "CEP192", "ESPL1", "FBLN7", "HMMR", "KIF18A",
    "KIF2C", "UBE2S", "GTSE1", "APOE", "CTNNA2", "EGF", "DOC2B", "PRSS21",
    "MIR100HG", "FGF19", "AP000439.2")
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = labelMarkers, transpose = TRUE)
svg("/storage/singlecell/zz4/fetal_bash/results/ArchR/PRPC/Plots/heatmapGS.svg")
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# addMotifAnnotations
PRPC <- addMotifAnnotations(ArchRProj = PRPC, motifSet = "cisbp", name = "Motif")
motifPositions <- getPositions(PRPC)
motifs <- c("ZFX")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
    value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
PRPC <- addGroupCoverages(ArchRProj = PRPC, groupBy = "Time_meta")
seFoot <- getFootprints(ArchRProj = PRPC, positions = motifPositions[markerMotifs],
    groupBy = "Time_meta")
plotFootprints(seFoot = seFoot, ArchRProj = PRPC, normMethod = "Divide",
    plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)

motifs <- labelMarkers
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions),
    value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
seFoot <- getFootprints(ArchRProj = PRPC, positions = motifPositions[markerMotifs],
    groupBy = "Time_meta")
plotFootprints(seFoot = seFoot, ArchRProj = PRPC, normMethod = "Divide",
    plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)

# getTrajectory()
trajMM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "MotifMatrix",
    log2Norm = FALSE)
mat <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),varCutOff = 0.7,returnMatrix = T)
subset_gene_names <- str_extract(rownames(mat), "(?<=:).*?(?=_)")
known_tfs <- c("AHR", "ASCL1", "ATOH7", "BARHL2", "BHLHE22", "CRX", "CTCF", 
"DLX1", "DLX2", "EBF1", "EBF2", "EBF3", "EGR1", "EOMES", "ESRRB", "FEZF1", "FEZF2", "FOXC1", "FOXN3", "FOXN4", "FOXP1", 
"HES1", "HES4", "ID1", "ID3", "IKZF2", "IKZF3", "INSM1", "INSM2", 
"IRX5", "IRX6", "ISL1", "LHX2", "LHX4", "LHX9", "MEF2C", "MYT1L", 
"NEUROD1", "NEUROD2", "NEUROD6", "NEUROG2", "NFIA", "NFIB", "NFIC", 
"NFIX", "NR2E1", "NR2E3", "NR4A2", "NRF1", "NRL", "OLIG2", "ONECUT1", 
"ONECUT2", "OTX2", "PAX6", "POU2F1", "POU2F2", "POU3F1", "POU4F2", 
"PRDM1", "PRDM13", "PRDM8", "PROX1", "PTF1A", "REST", "RORB", 
"SALL1", "SALL3", "SAMD11", "SMAD2", "SMAD3", "SMAD7", "SOX11", 
"SOX4", "SOX8", "SOX9", "TBR1", "TBX2", "TBX5", "TFAP2A", "TFAP2B", 
"TGIF1", "TGIF2", "THRB")
labels <- rownames(mat)[subset_gene_names %in% known_tfs]
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),varCutOff = 0.7,labelMarkers = labels,labelTop=0)
plotPDF(p1, name = "MotifMatrix.pdf", ArchRProj = PRPC,
    addDOC = FALSE, width = 6, height = 8)

trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGSM <- trajGSM[rowData(trajGSM)$name %in% subset_strings,]
mat <- plotTrajectoryHeatmap(trajGSM,varCutOff = 0.05,returnMatrix = T)
gs <- rownames(mat)

trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajGSM <- trajGSM[rowData(trajGSM)$name %in% subset_strings,]
mat <- plotTrajectoryHeatmap(trajGSM,varCutOff = 0.1,returnMatrix = T)
ge <- rownames(mat)
common_gene <- intersect(gs,ge)

trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
trajGSM <- trajGSM[rownames(trajGSM) %in% common_gene,]
mat <- plotTrajectoryHeatmap(trajGSM,varCutOff = 0.0,returnMatrix = T)
subset_str <- sub(".*:", "", rownames(mat))
labels <- rownames(mat)[subset_str %in% known_tfs]
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "blueYellow"),varCutOff = 0.0,labelMarkers = labels,labelTop=0)
plotPDF(p2, name = "GeneExpressionMatrix.pdf", ArchRProj = PRPC,
    addDOC = FALSE, width = 6, height = 8)

trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGSM <- trajGSM[rownames(trajGSM) %in% common_gene,]
subset_str <- sub(".*:", "", rownames(mat))
labels <- rownames(mat)[subset_str %in% known_tfs]
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"),varCutOff = 0.0,labelMarkers = labels,labelTop=0,rowOrder = rownames(mat))
plotPDF(p2, name = "GeneScoreMatrix.pdf", ArchRProj = PRPC,
    addDOC = FALSE, width = 6, height = 8)


trajGSM <- getTrajectory(ArchRProj = PRPC, name = "MyeloidU", useMatrix = "PeakMatrix", log2Norm = TRUE)
mat <- plotTrajectoryHeatmap(trajGSM,returnMatrix = T)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "greenBlue"),labelTop)
plotPDF(p2, name = "PeakMatrix.pdf", ArchRProj = PRPC,
    addDOC = FALSE, width = 6, height = 8)

PRPC$Days = as.numeric(PRPC$Days)
PRPC <- addTrajectory(ArchRProj = PRPC, name = "MyeloidU", groupBy = "Days",
    trajectory = trajectory, embedding = "UMAP", force = TRUE)
p <- plotTrajectory(PRPC, trajectory = "MyeloidU", colorBy = "cellColData", name = "Days")
plotPDF(p, name = "Plot-MyeloidU-Traj-Heatmaps.pdf", ArchRProj = PRPC,
    addDOC = FALSE, width = 6, height = 8)
