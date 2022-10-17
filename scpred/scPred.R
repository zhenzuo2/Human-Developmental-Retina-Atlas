args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_results_path <- args[2]
output_figures_path <- args[3]
sample_id <- args[4]
reference_file <- args[5]
meta_file <- args[6]

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
seurat_object@meta.data$Region <- meta[sample_id, "Region"]
seurat_object@meta.data$Days <- meta[sample_id, "Days"]

svg(paste(output_figures_path, sample_id, "_scpred_prediction.svg", sep = ""))
DimPlot(object = seurat_object, group.by = "scpred_prediction", label = TRUE)
dev.off()

saveRDS(seurat_object, paste(output_results_path, sample_id, "_scPred.rds",
    sep = ""))