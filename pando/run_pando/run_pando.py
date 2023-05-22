args <- commandArgs(trailingOnly = TRUE)
# Running Pando for each cell type
library(Seurat)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)
library(foreach)
library(doParallel)
registerDoParallel(8)
set.seed(0)

input_rna_file = args[1]
input_atac_file = args[2]
meta_file = args[3]
cell_type = args[4]
feature_selection = args[5]
n_features = as.numeric(args[6])
n_cells = as.numeric(args[7])
output_dir = args[8]
parallel = args[9]
mode = args[10]

# Check if the results already exist or not
dir.create(output_dir, showWarnings = FALSE)
if (file.exists(paste(output_dir, cell_type, "_modules_meta_",
    mode, "_feature_selection_", feature_selection, ".csv", sep = ""))){
    quit(save="no")
}

rna <- readRDS(input_rna_file)
meta <- read.csv(meta_file)
if (mode == "Macula") {
    meta <- meta[meta$Region == "Macula", ]
}
if (mode == "Peripheral") {
    meta <- meta[meta$Region == "Peripheral", ]
}
cells <- meta$X
rna <- subset(rna, cells = cells)
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna)
rna <- ScaleData(rna)
rna <- FindVariableFeatures(rna, nfeatures = n_features)
if (as.logical(feature_selection)) {
    rna <- subset(rna, features = VariableFeatures(rna))
}

seurat_object <- readRDS(input_atac_file)
common_cells <- intersect(colnames(seurat_object), colnames(rna))
seurat_object <- subset(seurat_object, cells = common_cells)
rna <- subset(rna, cells = common_cells)
seurat_object[["RNA"]] <- rna@assays$RNA
seurat_object@meta.data <- cbind(seurat_object@meta.data, rna@meta.data[colnames(seurat_object),
    ])

print(seurat_object)

if (length(colnames(seurat_object)) > n_cells) {
    cells = sample(colnames(seurat_object), n_cells)
    seurat_object <- subset(seurat_object, cells = cells)
}

seurat_object <- initiate_grn(seurat_object)
data(motifs)
seurat_object <- find_motifs(seurat_object, pfm = motifs, genome = BSgenome.Hsapiens.UCSC.hg38)

seurat_object <- infer_grn(seurat_object, peak_to_gene_method = "Signac",
    method = "glm", parallel = as.logical(parallel))
saveRDS(seurat_object, paste(output_dir, cell_type, "_seurat_object_",
    mode, "_feature_selection_", feature_selection, ".rds", sep = ""))
seurat_object <- find_modules(seurat_object)
modules <- NetworkModules(seurat_object)
write.csv(modules@meta, paste(output_dir, cell_type, "_modules_meta_",
    mode, "_feature_selection_", feature_selection, ".csv", sep = ""))