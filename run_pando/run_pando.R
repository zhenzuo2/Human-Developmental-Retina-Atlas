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
n_features = as.numeric(args[5])
output_dir = args[6]
parallel = args[7]
mode = args[8]

print("input_rna_file")
print(input_rna_file)

print("input_atac_file")
print(input_atac_file)

print("meta_file")
print(meta_file)

print("cell_type")
print(cell_type)

print("n_features")
print(n_features)

print("output_dir")
print(output_dir)

print("mode")
print(mode)

dir.create(output_dir, showWarnings = FALSE)

rna <- readRDS(input_rna_file)
meta <- read.csv(meta_file)
if (mode == "Macula"){
    meta<- meta[meta$Region == "Macula",]
}
if (mode == "Peripheral"){
    meta<- meta[meta$Region == "Peripheral",]
}
cells <- meta$X
rna <- subset(rna, cells = cells)
DefaultAssay(rna)<- "RNA"
rna <- NormalizeData(rna)
rna <- ScaleData(rna)
rna <- FindVariableFeatures(rna, nfeatures = n_features)
rna <- subset(rna, features = VariableFeatures(rna))

seurat_object <- readRDS(input_atac_file)
common_cells <- intersect(colnames(seurat_object), colnames(rna))
seurat_object <- subset(seurat_object, cells = common_cells)
rna <- subset(rna, cells = common_cells)
seurat_object[["RNA"]] <- rna@assays$RNA
seurat_object@meta.data<-cbind(seurat_object@meta.data, rna@meta.data[colnames(seurat_object),])

print(seurat_object)

if (length(colnames(seurat_object))>20000){
    cells = sample(colnames(seurat_object), 20000)
    seurat_object <- subset(seurat_object, cells = cells)
}

seurat_object <- initiate_grn(seurat_object)
data(motifs)
seurat_object <- find_motifs(seurat_object, pfm = motifs, genome = BSgenome.Hsapiens.UCSC.hg38)

seurat_object <- infer_grn(seurat_object, peak_to_gene_method = "Signac",
    method = "glm", parallel = as.logical(parallel))
saveRDS(seurat_object, paste(output_dir, cell_type, "_seurat_object_", mode, ".rds", sep = ""))
seurat_object <- find_modules(seurat_object)
modules <- NetworkModules(seurat_object)
write.csv(modules@meta, paste(output_dir, cell_type, "_modules_meta_", mode, ".csv", sep = ""))