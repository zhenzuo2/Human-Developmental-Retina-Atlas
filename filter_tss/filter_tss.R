args <- commandArgs(trailingOnly = TRUE)
library(Signac)
library(Seurat)

set.seed(0)
input_atac_file = args[1]
output_file = args[2]

seurat_object <- readRDS(input_atac_file)
seurat_object <- subset(x = seurat_object, subset = nCount_peaks > 1000 &
    nCount_peaks < 1e+05 & TSS.enrichment > 3)

saveRDS(seurat_object, output_file)