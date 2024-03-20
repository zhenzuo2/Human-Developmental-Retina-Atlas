library(Seurat)
library(Signac)

rna <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_rna/merged_rna.rds")
atac <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_atac/atac.rds")

common_cells <- intersect(colnames(rna), colnames(atac))

rna <- subset(rna, cells = common_cells)
atac <- subset(atac, cells = common_cells)
atac[["RNA"]] <- rna@assays$RNA

meta <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv")
rownames(meta) <- meta$X
meta <- meta[common_cells,]
atac@meta.data<-cbind(atac@meta.data, meta[colnames(atac),])

saveRDS(atac, "/storage/chentemp/zz4/adult_dev_compare/results/pando_merged_seurat_object/pando_merged_seurat_object.rds")