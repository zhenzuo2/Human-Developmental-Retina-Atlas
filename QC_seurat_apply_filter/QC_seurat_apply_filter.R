args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_results_path <- args[2]
output_figures_path <- args[3]
sample_id <- args[4]

library(Seurat)
seurat_object <- readRDS(input_file)

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

svg(paste(output_figures_path, sample_id, "_after_filtered_seurat_QC.svg", sep = ""))
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3, log = TRUE)
dev.off()

saveRDS(seurat_object, paste(output_results_path, sample_id, "_nFeature_RNA_500_5000_MT_5_fitered.rds", sep = ""))