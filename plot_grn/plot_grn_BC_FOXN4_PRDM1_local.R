# Change gene name to plot more genes
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
setwd("/Users/zz4/Desktop/")

seurat_object <- readRDS("BC_Pando/BC_seurat_object_All.rds")
# Read pseudotime object
cds <- readRDS("BC_monocle3_DE_analysis/monocle3.rds")
# Extract normalized matrix
normalized_matrix <- normalized_counts(cds)
# Extract pseudotime
time <- pseudotime(cds, reduction_method = "UMAP")

seurat_object <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
                              min_genes_per_module = 1, rsq_thresh = 0.05)
seurat_object <- get_network_graph(seurat_object, graph_name = "full_graph")
seurat_object <- get_tf_network(seurat_object, tf = "PRDM1", graph = "full_graph")
grn_plot <- plot_tf_network(seurat_object, tf = "PRDM1")

grn_plot$layers[[3]] <- NULL

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
grn_plot <- grn_plot + geom_point(aes(x = x, y = y, size = SUM, color = TIME),
                                  shape = shape) + scale_color_gradient(low = "blue", high = "red") +
  scale_size(10)
grn_plot <- grn_plot + ggrepel::geom_label_repel(data = dat, aes(x, y,
                                                                 label = labels, color = TIME), box.padding = 1, point.padding = 0,
                                                 max.overlaps = Inf, segment.color = "grey50", min.segment.length = 0)
ggsave(filename = paste("PRDM1_grn_plot.svg", sep = ""), plot = grn_plot,
       scale = 1, width = 10, height = 10)
