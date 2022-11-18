args <- commandArgs(trailingOnly = TRUE)
input_meta_path <- args[1]
output_results_path <- args[2]

library(dplyr)
library(Seurat)
meta <- data.frame()

input_meta_file = read.csv(input_meta_path)
input_file <- input_meta_file$Samples
for (file in input_file) {
    seurat_object <- readRDS(file)
    meta <- bind_rows(meta, seurat_object@meta.data)
}

meta = meta[, colSums(is.na(meta)) == 0]
meta <- meta[meta$scpred_prediction %in% c("AC", "RGC", "MG", "Cone", "BC", "HC","Rod"),]
meta <- meta[meta$scpred_max>0.99,]

write.csv(meta, paste(output_results_path, "meta.csv", sep = ""))