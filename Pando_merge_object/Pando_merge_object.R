args <- commandArgs(trailingOnly = TRUE)
input_RNA_path <- args[1]
input_ATAC_path <- args[2]
output_file <- args[3]

library(Seurat)
library(Pando)
set.seed(0)

rna <- readRDS(input_RNA_path)
DefaultAssay(rna)<- "RNA"
rna <- NormalizeData(rna)
rna <- ScaleData(rna)
seurat_object <- readRDS(input_ATAC_path)
common_cells <- intersect(colnames(seurat_object), colnames(rna))
seurat_object <- subset(seurat_object, cells = common_cells)
rna <- subset(rna, cells = common_cells)
seurat_object[["RNA"]] <- rna@assays$RNA
seurat_object@meta.data<-cbind(seurat_object@meta.data, rna@meta.data[colnames(seurat_object),])

saveRDS(seurat_object, output_file)
