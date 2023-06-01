args <- commandArgs(trailingOnly = TRUE)
name = args[1]
library(monocle3)
library(Seurat)
library(dplyr)

cells <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class.csv")
output_dir = "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/"

set.seed(0)
counts <- data.table::fread(paste("/storage/singlecell/zz4/fetal_snakemake/results/count_matrx_with_adult/",
    name, ".csv", sep = ""))
rownames(counts) <- counts$V1
counts$V1 <- NULL
temp <- t(counts)
colnames(temp) <- rownames(counts)
seurat_object <- Seurat::CreateSeuratObject(temp)
seurat_object@meta.data$sampleid <- as.character(colnames(seurat_object) %in% cells$X)
expression_matrix <- seurat_object@assays$RNA@counts

cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(row.names = rownames(expression_matrix))
gene_annotation$gene_short_name <- rownames(expression_matrix)
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata,
    gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment 
cds <- align_cds(cds, alignment_group = 'sampleid')

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds, use_partition = FALSE)

## Step 6: Order cells
root_cells = colnames(seurat_object)[!(colnames(seurat_object) %in% cells$X)]

cds <- order_cells(cds, root_cells = root_cells)

cds

saveRDS(cds, paste(output_dir, name, "_wadult_monocle3.rds", sep = ""))
write.table(pseudotime(cds, reduction_method = "UMAP"), paste(output_dir,
    name, "_wadult_monocle3_pseudotime.csv", sep = ""))