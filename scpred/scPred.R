args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
input_file <- args[2]
reference_file <- args[3]
meta_file <- args[4]
output_results_path <- args[5]
output_figures_path <- args[6]

suppressMessages(library(Seurat))
suppressMessages(library(scPred))

dir.create(output_figures_path, showWarnings = F)
dir.create(output_results_path, showWarnings = F)
set.seed(0)

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