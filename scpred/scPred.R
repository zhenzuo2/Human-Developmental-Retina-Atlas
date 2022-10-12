args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_results_path <- args[2]
sample_id <- args[3]
reference_file <- args[4]
meta_file <- args[5]

suppressMessages(library(Seurat))
suppressMessages(library(scPred))

seurat_object <- readRDS(input_file)
reference <- readRDS(reference_file)

seurat_object <- NormalizeData(seurat_object)
seurat_object <- scPredict(new = seurat_object, reference = reference,
    max.iter.harmony = 100)
DefaultAssay(seurat_object) <- "RNA"

meta <- read.csv(meta_file)
rownames(meta) <- meta$Samples

seurat_object@meta.data$Time <- meta[sample_id, "Time"]
seurat_object@meta.data$Reion <- meta[sample_id, "Reion"]
seurat_object@meta.data$Days <- meta[sample_id, "Days"]

saveRDS(seurat_object, paste(output_results_path, sample_id, "_scPred.rds", sep = ""))