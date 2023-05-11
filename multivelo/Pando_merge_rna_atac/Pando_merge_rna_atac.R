input_RNA_path="/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds"
input_ATAC_path="/storage/singlecell/zz4/fetal_snakemake/results/merged_atac/atac.rds"
output_file="/storage/singlecell/zz4/fetal_snakemake/results/Pando_merged/seurat_object.rds"

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
