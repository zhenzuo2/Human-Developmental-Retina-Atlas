seurat_pando_object <- "/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_All_feature_selection_FALSE.rds"
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(monocle3)

set.seed(0)
# Read pando result rds object
seurat_object <- readRDS(seurat_pando_object)
# Read pseudotime object
cds <- readRDS(time_cds)
# Extract normalized matrix
normalized_matrix <- normalized_counts(cds)
# Extract pseudotime
time <- pseudotime(cds, reduction_method = "UMAP")
# Find and filter modles
seurat_object <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
    min_genes_per_module = 1, rsq_thresh = 0.05)
modules <- NetworkModules(seurat_object)
modules
# plot quality control figures
gof <- plot_gof(seurat_object, point_size = 3)
ggsave(filename = paste(output_dir, label, "_gof.svg", sep = ""), plot = gof, scale = 1,
    width = 10, height = 10)

module_metrics <- plot_module_metrics(seurat_object)
ggsave(filename = paste(output_dir, label, "_module_metrics.svg", sep = ""), plot = module_metrics,
    scale = 1, width = 10, height = 10)

# Get network graph
seurat_object <- get_network_graph(seurat_object, graph_name = "umap_graph",
    umap_method = "weighted")
# Save the initial plot
grn_plot <- plot_network_graph(seurat_object, graph = "umap_graph", label_nodes = FALSE)
grn_plot$layers[[2]] <- NULL

# Find all TFs in the GRN
TFS <- unique(seurat_object@grn@networks$glm_network@modules@meta$tf)

dat <- grn_plot$data
dat$TIME = 0
dat$SUM = 0
rownames(dat) = dat$name
TFS <- intersect(TFS, dat$name)

for (gene in dat$name) {
    dat[gene, "SUM"] = sum(normalized_matrix[gene, ])
    dat[gene, "TIME"] = sum(time * normalized_matrix[gene, ])/sum(normalized_matrix[gene,
        ])
}
dat$TIME <- scales::squish(dat$TIME, quantile(dat$TIME, c(0.01, 0.99)))

dat$is_TF <- dat$name %in% TFS
dat$labels <- ""
dat[dat$name[dat$is_TF], "labels"] <- dat$name[dat$is_TF]
dat$shape <- ifelse(dat$is_TF, 19, 1)
attach(dat)
grn_plot <- grn_plot + geom_point(aes(x = UMAP_1, y = UMAP_2, size = SUM,
    color = TIME), shape = shape) + scale_color_gradient(low = "blue",
    high = "red") + scale_size(10)
grn_plot <- grn_plot + ggrepel::geom_label_repel(data = dat, aes(UMAP_1,
    UMAP_2, label = labels, color = TIME), box.padding = 1, point.padding = 0,
    max.overlaps = Inf, segment.color = "grey50", min.segment.length = 0)
ggsave(filename = paste(output_dir, label, "_grn_plot.svg", sep = ""), plot = grn_plot,
    scale = 1, width = 10, height = 10)
write.csv(dat,paste(output_dir, label, "_grn_plot.csv", sep = ""))