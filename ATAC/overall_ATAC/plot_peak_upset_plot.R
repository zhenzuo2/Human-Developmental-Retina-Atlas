output_results_path = "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-ALL")

p2g <- getPeak2GeneLinks(ArchRProj = proj2, corCutOff = 0.45, resolution = 1,
    returnLoops = FALSE)

p2g <- as.data.frame(metadata(p2g)[[1]][p2g$idxATAC, ])

p2g <- p2g[!duplicated(p2g), ]

p2g$ID <- paste(p2g$seqnames, p2g$start, p2g$end, sep = "_")

allpeaks <- data.frame(proj2@peakSet)
allpeaks$ID <- paste(allpeaks$seqnames, allpeaks$start, allpeaks$end, sep = "_")

markersPeaks <- getMarkerFeatures(ArchRProj = proj2, useMatrix = "PeakMatrix",
    groupBy = "majorclass", bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
res <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(res) <- c("ID", "celltypes")
for (celltypes in c("AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC",
    "Rod")) {
    ID <- paste(markerList[[celltypes]]$seqnames, "_", markerList[[celltypes]]$start,
        "_", markerList[[celltypes]]$end, sep = "")
    res <- rbind(res, data.frame(ID, celltypes))
}

AC <- res[res$celltypes == "AC", "ID"]
BC <- res[res$celltypes == "BC", "ID"]
Cone <- res[res$celltypes == "Cone", "ID"]
HC <- res[res$celltypes == "HC", "ID"]
MG <- res[res$celltypes == "MG", "ID"]
NRPC <- res[res$celltypes == "NRPC", "ID"]
PRPC <- res[res$celltypes == "PRPC", "ID"]
RGC <- res[res$celltypes == "RGC", "ID"]
Rod <- res[res$celltypes == "Rod", "ID"]
# Load the required library
library(UpSetR)

# Create your sets of strings
All_Peaks <- allpeaks$ID
Linked_Peaks <- p2g$ID
DAR <- res$ID

# Combine the sets into a list
sets_list <- list(All_Peaks = All_Peaks, Linked_Peaks = Linked_Peaks, AC = AC,
    BC = BC, Cone = Cone, HC = HC, MG = MG, NRPC = NRPC, PRPC = PRPC, RGC = RGC,
    Rod = Rod)
intersections = list(list("Linked_Peaks"), list("Linked_Peaks", "AC"), list("Linked_Peaks", "BC"),
    list("Linked_Peaks", "Cone"), list("Linked_Peaks", "HC"), list("Linked_Peaks",
        "MG"), list("Linked_Peaks", "NRPC"), list("Linked_Peaks", "PRPC"),
    list("Linked_Peaks", "RGC"), list("Linked_Peaks", "Rod"))
# Create the upset plot
upset_data <- fromList(sets_list)
upset_plot <- upset(upset_data, order.by = "freq", main.bar.color = "#56B4E9",
    nsets = 11, intersections = intersections,text.scale = 2)

# Display the plot
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp6.svg", width = 10,
    height = 10)
print(upset_plot)
dev.off()
