args <- commandArgs(trailingOnly = TRUE)
library(monocle3)
library(Seurat)
library(dplyr)

input_file = args[1]
meta_file = args[2]
output_dir = args[3]
mode = args[4]

set.seed(0)
seurat_object <- readRDS(input_file)
meta <- read.csv(meta_file)
meta <- meta[meta$Region %in% c("Macula", "Peripheral"), ]
rownames(meta) <- meta$X
cells <- meta$X
common_cells <- intersect(cells, colnames(seurat_object))
seurat_object <- subset(seurat_object, cells = common_cells)
meta <- meta[common_cells,]

seurat_object@meta.data$sampleid <- meta[rownames(seurat_object@meta.data),
    "sampleid"]
seurat_object@meta.data$Days <- meta[rownames(seurat_object@meta.data),
    "Days"]
seurat_object@meta.data$Region <- meta[rownames(seurat_object@meta.data),
    "Region"]

expression_matrix <- seurat_object@assays$RNA@counts

cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(row.names = rownames(expression_matrix))
gene_annotation$gene_short_name <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata,
    gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

region_days_models <- fit_models(cds, model_formula_str = "~ Region + Days",
    verbose = TRUE, expression_family = "negbinomial")
fit_coefs_region_days_models <- coefficient_table(region_days_models)

days_models <- fit_models(cds, model_formula_str = "~ Days", expression_family = "negbinomial",
    verbose = TRUE)
fit_coefs_days_models <- coefficient_table(days_models)

region_models <- fit_models(cds, model_formula_str = "~ Region", expression_family = "negbinomial",
    verbose = TRUE)
fit_coefs_region_models <- coefficient_table(region_models)

compare_mod_region <- compare_models(region_days_models, days_models)

compare_mod_days <- compare_models(region_days_models, region_models)

write.table(fit_coefs_region_days_models[, c(-4, -5)], paste(output_dir,
    mode, "_fit_coefs_region_days_models.csv", sep = ""), quote = FALSE,
    sep = ",", row.names = FALSE)
write.table(fit_coefs_days_models[, c(-4, -5)], paste(output_dir, mode,
    "_fit_coefs_days_models.csv", sep = ""), quote = FALSE, sep = ",",
    row.names = FALSE)
write.table(fit_coefs_region_models[, c(-4, -5)], paste(output_dir, mode,
    "_fit_coefs_region_models.csv", sep = ""), quote = FALSE, sep = ",",
    row.names = FALSE)

write.table(compare_mod_region, paste(output_dir, mode, "_compare_mod_region.csv",
    sep = ""), quote = FALSE, sep = ",", row.names = FALSE)
write.table(compare_mod_days, paste(output_dir, mode, "_compare_mod_days.csv",
    sep = ""), quote = FALSE, sep = ",", row.names = FALSE)