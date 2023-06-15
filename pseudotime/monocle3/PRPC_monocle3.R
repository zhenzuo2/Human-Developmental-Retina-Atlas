args <- commandArgs(trailingOnly = TRUE)
library(monocle3)
library(Seurat)
library(dplyr)

input_file = "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds"
meta_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
output_dir = "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/"
name = "PRPC"

set.seed(0)
seurat_object <- readRDS(input_file)
meta <- read.csv(meta_file)
meta <- meta[meta$Region %in% c("Macula", "Peripheral"), ]
meta <- meta[meta$majorclass %in% c("PRPC"), ]
rownames(meta) <- meta$X
cells <- meta$X
common_cells <- intersect(cells, colnames(seurat_object))
seurat_object <- subset(seurat_object, cells = common_cells)
meta <- meta[common_cells, ]

seurat_object@meta.data$sampleid <- meta[rownames(seurat_object@meta.data),
    "sampleid"]

seurat_object <- NormalizeData(seurat_object)

expression_matrix <- seurat_object@assays$RNA@counts

cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(row.names = rownames(expression_matrix))
gene_annotation$gene_short_name <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata,
    gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "sampleid")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = FALSE)

## Step 6: Order cells
root_cells <- colnames(cds)[cds$sampleid == "Multiome_10w_NR"]
if (length(root_cells) == 0) {
    root_cells <- colnames(cds)[cds$Days == min(unique(cds$Days))]
}

cds <- order_cells(cds, root_cells = root_cells)

cds

saveRDS(cds, paste(output_dir, name, "_monocle3.rds", sep = ""))
write.table(pseudotime(cds, reduction_method = "UMAP"), paste(output_dir,
    name, "_monocle3_pseudotime.csv", sep = ""))