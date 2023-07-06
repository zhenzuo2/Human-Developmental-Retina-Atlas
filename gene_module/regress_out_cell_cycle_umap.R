library(Seurat)
library(Matrix)

seurat_object <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/")
seurat_object@meta.data$sampleid = sub("_[^_]*$", "", rownames(seurat_object@meta.data))
# split the dataset into a list of two seurat objects (stim and CTRL)
samples <- names(table(seurat_object@meta.data$sampleid)[table(seurat_object@meta.data$sampleid)>100])
cells <- rownames(seurat_object@meta.data[seurat_object@meta.data$sampleid %in% samples,])
seurat_object <- subset(seurat_object,cells = cells)
seurat_object.list <- SplitObject(seurat_object, split.by = "sampleid")

# normalize and identify variable features for each dataset independently
seurat_object.list <- lapply(X = seurat_object.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat_object.list)

seurat_object.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, anchor.features = features)

# this command creates an 'integrated' data assay
seurat_object.combined <- IntegrateData(anchorset = seurat_object.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat_object.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_object.combined <- ScaleData(seurat_object.combined, verbose = FALSE)
seurat_object.combined <- RunPCA(seurat_object.combined, npcs = 30, verbose = FALSE)
seurat_object.combined <- RunUMAP(seurat_object.combined, reduction = "pca", dims = 1:30)
seurat_object.combined <- FindNeighbors(seurat_object.combined, reduction = "pca", dims = 1:30)
seurat_object.combined <- FindClusters(seurat_object.combined, resolution = 0.5)