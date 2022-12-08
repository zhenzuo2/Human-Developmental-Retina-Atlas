args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
meta_file <- args[2]
sample_name <- args[3]
output_path <- args[4]

dir.create(output_path, showWarnings = FALSE)


# MultiVelo Seurat WNN Demo
# The procedure mostly follows Seurat tutorial: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# Note that we do not claim these preprocessing steps to be the best, as there aren't any. Feel free to make any changes you deem necessary.
# Please use libblas 3.9.1 and liblapack 3.9.1 for reproducing the 10X mouse seurat_object demo, or use supplied WNN files on GitHub.

library(Seurat)
library(Signac)
# read in expression and accessbility data
seurat_object <- readRDS(input_file)

# subset for the same cells in the jointly filtered anndata object
barcodes <- read.csv(meta_file)

filtered_cells <- read.csv("/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/filtered_cells.csv")
seurat_object <- subset(seurat_object, cells = intersect(barcodes$X,filtered_cells$X))
seurat_object

# preprocess RNA
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object) # not scaled for consistency with scVelo (optionally, use SCTransform)
seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional

# preprocess ATAC
DefaultAssay(seurat_object) <- "peaks"
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object)
seurat_object <- RunSVD(seurat_object)
seurat_object <- RunUMAP(seurat_object, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional

# find weighted nearest neighbors
seurat_object <- FindMultiModalNeighbors(seurat_object, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 50)
seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional

# extract neighborhood graph
nn_idx <- seurat_object@neighbors$weighted.nn@nn.idx
nn_dist <- seurat_object@neighbors$weighted.nn@nn.dist
nn_cells <- seurat_object@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx, paste(output_path, sample_name, "_nn_idx.txt",sep=""), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, paste(output_path, sample_name, "_nn_dist.txt",sep=""), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, paste(output_path, sample_name, "_nn_cells.txt",sep=""), sep = ',', row.names = F, col.names = F, quote = F)
saveRDS(seurat_object, paste(output_path, sample_name, "_seurat_object.rds",sep=""))