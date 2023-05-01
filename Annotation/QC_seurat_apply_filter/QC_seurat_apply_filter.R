args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
input_file <- args[2]
output_results_path <- args[3]
output_figures_path <- args[4]

library(Seurat)
set.seed(0)
dir.create(output_figures_path, showWarnings = F)
dir.create(output_results_path, showWarnings = F)

seurat_object <- readRDS(input_file)

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10)

svg(paste(output_figures_path, sample_id, "_after_filtered_seurat_QC.svg", sep = ""))
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3, log = TRUE)
dev.off()

saveRDS(seurat_object, paste(output_results_path, sample_id, ".rds", sep = ""))