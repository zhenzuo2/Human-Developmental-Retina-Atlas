library(Seurat)
input_rna_file = "/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
seurat_object <- readRDS(input_rna_file)

cells <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC.csv")
seurat_object <- subset(seurat_object, cells = cells$X)

saveRDS(seurat_object,"/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna_PRPC.rds")