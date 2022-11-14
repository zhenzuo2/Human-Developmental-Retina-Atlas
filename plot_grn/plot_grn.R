library(Pando)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

seurat_object <- readRDS(input_file)
modules <- NetworkModules(seurat_object)
seurat_object <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
    min_genes_per_module = 1, rsq_thresh = 0.05)
plot_gof(seurat_object, point_size = 3)

svg(paste(output_dir, "Pando_QC.svg"))
plot_module_metrics(seurat_object)
dev.off()

seurat_object <- get_network_graph(seurat_object, graph_name = "umap_graph")

svg(paste(output_dir, "Pando_TF_TG.svg"))
plot_network_graph(seurat_object, graph = "umap_graph")
dev.off()

modules <- NetworkModules(seurat_object)
seurat_object <- get_network_graph(seurat_object, graph_name = "sub_graph",
    umap_method = "corr", features = unique(modules@meta$tf))

svg(paste(output_dir, "Pando_TF_TG.svg"))
plot_network_graph(seurat_object, graph = "sub_graph")
dev.off()