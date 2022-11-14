args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
input_file <- args[2]
output_results_path <- args[3]
output_figures_path <- args[4]


# load libraries
library(DoubletFinder)
library(Seurat)
seurat_object <- readRDS(input_file)
# https://rpubs.com/kenneditodd/doublet_finder_example

## Pre-process Seurat object (standard)
## --------------------------------------------------------------------------------------
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",
    nfeatures = 2000)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)

# Find significant PCs
stdv <- seurat_object[["pca"]]@stdev
sum.stdv <- sum(seurat_object[["pca"]]@stdev)
percent.stdv <- (stdv/sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) >
    0.1), decreasing = T)[1] + 1
min.pc <- min(co1, co2)

print("min.pc")
print(min.pc)

# finish pre-processing
seurat_object <- RunUMAP(seurat_object, dims = 1:min.pc)
seurat_object <- FindNeighbors(object = seurat_object, dims = 1:min.pc)
seurat_object <- FindClusters(object = seurat_object, resolution = 0.1)

# pK identification (no ground-truth)
sweep.list <- paramSweep_v3(seurat_object, PCs = 1:min.pc)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)

# Optimal pK is the max of the bomodality coefficent (BCmvn)
# distribution
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
print("optimal.pk")
print(optimal.pk)
print("optimal.pk")
print(optimal.pk)

## Homotypic doublet proportion estimate
annotations <- seurat_object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp.poi <- round(optimal.pk * nrow(seurat_object@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
seurat_object <- doubletFinder_v3(seu = seurat_object, PCs = 1:min.pc,
    pK = optimal.pk, nExp = nExp.poi.adj)
metadata <- seurat_object@meta.data
colnames(metadata)[ncol(metadata)] <- "doublet_finder"
seurat_object@meta.data <- metadata

svg(paste(output_figures_path, sample_id, "_log_seurat_QC.svg", sep = ""))
DimPlot(object = seurat_object,group.by = "doublet_finder")
dev.off()

# subset and save
seurat_object <- subset(seurat_object, doublet_finder == "Singlet")

saveRDS(seurat_object, paste(output_results_path, sample_id, ".rds", sep = ""))