# Take all args
args <- commandArgs(trailingOnly = TRUE)
library(Signac)
library(Seurat)

set.seed(0)
input_atac_file = args[1]
output_file = args[2]

seurat_object <- readRDS(input_atac_file)
seurat_object <- TSSEnrichment(seurat_object, fast = FALSE)

saveRDS(seurat_object, output_file)
