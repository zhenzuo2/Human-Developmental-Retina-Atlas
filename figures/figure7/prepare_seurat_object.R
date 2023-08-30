library(Seurat)
seurat_object <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cellbygene/adata.obs.csv")

seurat_object_ <- subset(seurat_object,cells = meta$X)

rownames(meta) <- meta$X

SaveH5Seurat(seurat_object_,"/storage/singlecell/zz4/fetal_snakemake/results/cellbygene/adata.h5Seurat")