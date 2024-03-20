output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
library(GenomicRanges)
set.seed(0)
setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj3")

# Define the ranges as a character vector
ranges <- c("chr14:56888304-56888662", "chr14:56964627-56965216", "chr14:56902629-56903105",
    "chr14:56899344-56901557", "chr14:56820328-56822794", "chr14:56807853-56808469",
    "chr14:56805457-56806309", "chr14:56805457-56807186", "chr14:56730099-56731208",
    "chr14:56742261-56743173", "chr14:56754441-56756016", "chr14:56735566-56735995",
    "chr14:56753241-56754311")

nams <- c("FM1", "FM3", "AN1", "AN2", "EELPOT", "ECR2", "CM", "VE", "O5",
    "O7", "O9", "O6", "O8")

# Create a GRanges object
gr <- GRanges(ranges)
nams <- nams[order(gr)]
# Print the GRanges object
print(gr)

markerGenes <- c("OTX2")
p <- plotBrowserTrack(ArchRProj = proj2, features = gr, groupBy = "majorclass",
    useGroups = c("BC", "Rod", "Cone", "NRPC", "AC", "MG", "HC", "PRPC",
        "RGC"), geneSymbol = markerGenes, loops = getPeak2GeneLinks(proj2),
    upstream = 1e+05, downstream = 2e+05, )
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf",
    ArchRProj = proj2, addDOC = FALSE, width = 8, height = 5)
