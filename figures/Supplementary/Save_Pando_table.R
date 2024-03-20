########################################################################
#AC
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggforce)

output_dir = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"
label = "AC"
seurat_pando_object = "/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_seurat_object_All_feature_selection_FALSE.rds"
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

tfs <- dat$name[dat$name %in% TFS]
temp <- data.frame(coef(seurat_object))
temp <- temp[temp$tf %in% tfs,]
write.csv(coef(seurat_object),"/storage/singlecell/zz4/fetal_snakemake/temp/AC.csv")

write.csv(temp,"/storage/singlecell/zz4/fetal_snakemake/temp/AC2.csv")

########################################################################
#PRPC
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggforce)

output_dir = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
label = "PRPC"
seurat_pando_object = "/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_All_feature_selection_FALSE.rds"
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

tfs <- dat$name[dat$name %in% TFS]
temp <- data.frame(coef(seurat_object))
temp <- temp[temp$tf %in% tfs,]
write.csv(coef(seurat_object),"/storage/singlecell/zz4/fetal_snakemake/temp/PRPC.csv")

write.csv(temp,"/storage/singlecell/zz4/fetal_snakemake/temp/PRPC2.csv")


########################################################################
#NRPC
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
tfs <- dat$name[dat$name %in% TFS]
temp <- data.frame(coef(seurat_object))
temp <- temp[temp$tf %in% tfs,]
write.csv(coef(seurat_object),"/storage/singlecell/zz4/fetal_snakemake/temp/BC.csv")

write.csv(temp,"/storage/singlecell/zz4/fetal_snakemake/temp/BC2.csv")
