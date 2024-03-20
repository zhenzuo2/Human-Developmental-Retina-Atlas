args <- commandArgs(trailingOnly = TRUE)

n_cells <- 20000
output_dir <- "/storage/chentemp/zz4/adult_dev_compare/results/pando_results/pando_"

# Running Pando for each cell type
library(Seurat)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
library(foreach)
library(doParallel)
registerDoParallel(8)
set.seed(0)

input_file <- args[1]
cell_type <- args[2]
seurat_object <- readRDS(input_file)

if (length(colnames(seurat_object)) > n_cells) {
    cells <- sample(colnames(seurat_object), n_cells)
    seurat_object <- subset(seurat_object, cells = cells)
}

seurat_object <- initiate_grn(seurat_object)
data(motifs)
seurat_object <- find_motifs(seurat_object, pfm = motifs, genome = BSgenome.Hsapiens.UCSC.hg38)
seurat_object <- infer_grn(seurat_object, peak_to_gene_method = "Signac",
    method = "glm", parallel = TRUE)
saveRDS(seurat_object, paste(output_dir, cell_type, ".rds", sep = ""))
