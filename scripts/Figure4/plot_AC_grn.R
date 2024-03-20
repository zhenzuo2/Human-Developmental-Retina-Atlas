output_dir = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/"
label = "AC"
dir.create(output_dir, showWarnings = FALSE)
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)

set.seed(0)
# Read pando result rds object
seurat_object <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/pando_results/pando_AC.rds")
# Extract normalized matrix
normalized_matrix <- seurat_object@assays$RNA@data
# Extract pseudotime
time <- seurat_object@meta.data$Days
names(time) <- rownames(seurat_object@meta.data)
# Find and filter modles
seurat_object <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
    min_genes_per_module = 1, rsq_thresh = 0.05)
modules <- NetworkModules(seurat_object)

modules
write.csv(modules@meta, paste(output_dir, label, "_grn_meta.csv", sep = ""))

# plot quality control figures
gof <- plot_gof(seurat_object, point_size = 3)
ggsave(filename = paste(output_dir, label, "_gof.tiff", sep = ""), plot = gof,
    scale = 1, width = 10, height = 10)

module_metrics <- plot_module_metrics(seurat_object)
ggsave(filename = paste(output_dir, label, "_module_metrics.tiff", sep = ""),
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
dat$TIME <- scales::squish(dat$TIME, quantile(dat$TIME, c(0.01, 0.99)))
dat$is_TF <- dat$name %in% TFS
dat$labels <- ""
dat[dat$name[dat$is_TF], "labels"] <- dat$name[dat$is_TF]
dat$shape <- ifelse(dat$is_TF, 19, 1)
dat$SUM = dat$SUM/max(dat$SUM)
dat$Nor_Exp = dat$SUM*10
write.csv(dat, paste(output_dir, label, "_grn_plot.csv", sep = ""))
##########################################################################################
p <- ggplot(dat, aes(x, y))
p <- p + geom_point(aes(x = UMAP_1, y = UMAP_2, color = TIME,size = Nor_Exp),
    shape = dat$shape) + scale_color_gradient(low = "blue", high = "red") +
    scale_size(10)
p <- p + ggrepel::geom_text_repel(data = dat, size = 8, aes(UMAP_1, UMAP_2,
    label = labels, color = TIME), box.padding = 1, point.padding = 0,
    max.overlaps = Inf, segment.color = "grey50", min.segment.length = 0)
p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) + theme(panel.background = element_blank())
p <- p +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
ggsave(filename = paste(output_dir, label, "_grn_plot.tiff", sep = ""),
    plot = p, scale = 1, width = 24, height = 18)

muo_data <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
    min_genes_per_module = 1, rsq_thresh = 0.05)
muo_data <- get_network_graph(muo_data, graph_name = "full_graph", umap_method = "none")

muo_data <- get_tf_network(muo_data, tf = "ONECUT1", graph = "full_graph")
p <- plot_tf_network(muo_data, tf = "ONECUT1")
ggsave(filename = paste(output_dir, label, "ONECUT1.tiff", sep = ""), plot = p,
    scale = 1, width = 10, height = 10)

k <- 2  # Number of clusters
kmeans_result <- kmeans(dat[, c("x", "y")], centers = k)

# Add cluster labels to the data frame
dat$cluster <- as.factor(kmeans_result$cluster)
dat <- dat[dat$labels!="",]
p <- ggplot(dat, aes(x, y))
p <- p + geom_point(aes(x = UMAP_1, y = UMAP_2, color = cluster))+scale_size(10)
ggsave(filename = paste(output_dir, label, "_grn_plot_cluster.tiff", sep = ""),
    plot = p, scale = 1, , width = 18, height = 18)
write.csv(dat,"/storage/chentemp/zz4/adult_dev_compare/temp/temp.csv")
