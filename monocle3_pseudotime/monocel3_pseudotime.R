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

seurat_object@meta.data$batch <- meta[rownames(seurat_object@meta.data),
    "batch"]
seurat_object@meta.data$scpred_prediction <- meta[rownames(seurat_object@meta.data),
    "scpred_prediction"]
seurat_object@meta.data$Time <- meta[rownames(seurat_object@meta.data),
    "Time"]
seurat_object@meta.data$Region <- meta[rownames(seurat_object@meta.data),
    "Region"]
seurat_object@meta.data$Days <- meta[rownames(seurat_object@meta.data),
    "Days"]
seurat_object@meta.data$leiden <- meta[rownames(seurat_object@meta.data),
    "leiden"]
seurat_object@meta.data$scpred_prediction_mode <- meta[rownames(seurat_object@meta.data),
    "scpred_prediction_mode"]
seurat_object@meta.data$subclass <- meta[rownames(seurat_object@meta.data),
    "subclass"]
seurat_object@meta.data$majorclass <- meta[rownames(seurat_object@meta.data),
    "majorclass"]

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",
    nfeatures = 10000)

expression_matrix <- seurat_object@assays$RNA@counts[VariableFeatures(seurat_object),
    ]
cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(row.names = rownames(expression_matrix))
gene_annotation$gene_short_name <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata,
    gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
if ((meta_file != "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_obs.csv") &
    (meta_file != "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_obs.csv")) {
    cds <- align_cds(cds, alignment_group = "batch")
}

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = FALSE)

## Step 6: Order cells
root_cells <- colnames(cds)[cds$batch == "Multiome_10w_NR"]
if (length(root_cells) == 0) {
    root_cells <- colnames(cds)[cds$Days == 70]
}

if (meta_file == "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_obs.csv") {
    root_cells <- colnames(cds)[cds$Days == 62]
}
print("Root_cells")
print(root_cells)
cds <- order_cells(cds, root_cells = root_cells)

saveRDS(cds, paste(output_dir, "monocle3.rds", sep = ""))
write.table(pseudotime(cds, reduction_method = "UMAP"), paste(output_dir,
    "monocle3_pseudotime.csv", sep = ""))
