# Take all args
args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
sample_id <- args[2]
min_cells <- args[3]
min_features <- args[4]
output_figures_path <- args[5]
output_results_path <- args[6]

# Read packages
library(dplyr)
library(Seurat)
library(patchwork)

# Read h5 files
counts <- Read10X_h5(paste0(input_path, sample_id, "/outs/filtered_feature_bc_matrix.h5",
    sep = ""))

# Create a Seurat object containing the RNA data
seurat_object <- CreateSeuratObject(counts = counts$`Gene Expression`,
    assay = "RNA", min.cells = min_cells, min.features = min_features)
DefaultAssay(seurat_object) <- "RNA"

# Rename cells with sample ID
seurat_object <- RenameCells(object = seurat_object, add.cell.id = sample_id)
seurat_object <- NormalizeData(seurat_object)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
DefaultAssay(seurat_object) <- "RNA"

# Plot QC images
svg(paste(output_figures_path, sample_id, "_seurat_QC.svg", sep = ""))
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3)
dev.off()

# Plot log QC images
svg(paste(output_figures_path, sample_id, "_log_seurat_QC.svg", sep = ""))
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3, log = TRUE)
dev.off()

# Save RDS
saveRDS(seurat_object, paste(output_results_path, sample_id, "_min_cell_",
    as.character(min_cells), "_min_features_", as.character(min_features),
    ".rds", sep = ""))