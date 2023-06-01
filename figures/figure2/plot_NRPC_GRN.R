library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggforce)

output_dir = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"
label = "NRPC"
seurat_pando_object = "/storage/singlecell/zz4/fetal_snakemake/results/pando_results/Annotated_NRPC_seurat_object.rds"
time_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv"
set.seed(0)
# Read pando result rds object
seurat_object <- readRDS(seurat_pando_object)
# Read pseudotime object
df <- read.csv(time_file)
rownames(df) <- df$X
# Extract normalized matrix
normalized_matrix <- seurat_object@assays$RNA@data
cells <- intersect(rownames(df), colnames(normalized_matrix))
normalized_matrix <- normalized_matrix[, cells]
# Extract pseudotime
time <- df[cells, "Days"]

# Find and filter modles
seurat_object <- find_modules(seurat_object)
modules <- NetworkModules(seurat_object)
modules
# plot quality control figures
gof <- plot_gof(seurat_object, point_size = 3)
ggsave(filename = paste(output_dir, label, "_gof.svg", sep = ""), plot = gof,
       scale = 1, width = 10, height = 10)

module_metrics <- plot_module_metrics(seurat_object)
ggsave(filename = paste(output_dir, label, "_module_metrics.svg", sep = ""),
       plot = module_metrics, scale = 1, width = 10, height = 10)

# Get network graph
seurat_object <- get_network_graph(seurat_object, graph_name = "umap_graph",
                                   umap_method = "weighted")
# Save the initial plot
grn_plot <- plot_network_graph(seurat_object, graph = "umap_graph", label_nodes = FALSE)

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

dat$is_TF <- dat$name %in% TFS
dat$labels <- ""
dat[dat$name[dat$is_TF], "labels"] <- dat$name[dat$is_TF]
dat$shape <- ifelse(dat$is_TF, 19, 1)
dat = dat[dat$labels != "", ]


k <- 3  # Number of clusters
kmeans_result <- kmeans(dat[, c("x", "y")], centers = k)

# Add cluster labels to the data frame
dat$cluster <- as.factor(kmeans_result$cluster)

grn_plot_ <- ggplot() + geom_point(data = dat, aes(x = UMAP_1, y = UMAP_2,
                                                   size = SUM, color = TIME), shape = dat$shape) + scale_color_gradient(low = "blue",
                                                                                                                        high = "red") + scale_size(5)
grn_plot_final <- grn_plot_ + ggrepel::geom_label_repel(data = dat, size = 10,
                                                        aes(UMAP_1, UMAP_2, label = labels, color = TIME), box.padding = 1,
                                                        point.padding = 0, label.size = 0, max.overlaps = Inf, segment.color = "grey50",
                                                        min.segment.length = 0, fill = NA) + theme(plot.background = element_blank(),
                                                                                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                   panel.border = element_blank()) + theme(axis.title.x = element_blank(),
                                                                                                                                           axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(),
                                                                                                                                           axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(panel.background = element_rect(fill = "transparent"),
                                                                                                                                                                                                                  plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(),
                                                                                                                                                                                                                  panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"),
                                                                                                                                                                                                                  legend.box.background = element_rect(fill = "transparent"))
ggsave(filename = paste(output_dir, label, "_grn_plot.svg", sep = ""),
       plot = grn_plot_final, scale = 1, width = 20, height = 15, bg = "transparent")

grn_plot_final2 <- ggplot(dat, aes(x = x, y = y)) + geom_point(aes(color = as.factor(cluster))) +
  geom_mark_circle(aes(color = as.factor(cluster)), expand = unit(0.5,
                                                                  "mm")) + theme(legend.position = "none")

ggsave(filename = paste(output_dir, label, "_grn_plot_cluster.svg", sep = ""),
       plot = grn_plot_final2, scale = 1, width = 20, height = 15)
write.csv(dat, paste(output_dir, label, "_grn_plot.csv", sep = ""))
