library(monocle3)
library(Seurat)
library(dplyr)
library(SCENT)
data(net13Jun12)
print(dim(net13Jun12.m))
library("org.Hs.eg.db")

input_file = "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds"
meta_file = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
output_dir = "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/"
name = "PRPC_MG"

set.seed(0)
seurat_object <- readRDS(input_file)
meta <- read.csv(meta_file)
meta <- meta[meta$Region %in% c("Macula", "Peripheral"), ]
meta <- meta[meta$majorclass %in% c("MG", "PRPC"), ]
rownames(meta) <- meta$X
cells <- meta$X
common_cells <- intersect(cells, colnames(seurat_object))
seurat_object <- subset(seurat_object, cells = common_cells)
meta <- meta[common_cells, ]

seurat_object <- NormalizeData(seurat_object)

counts <- seurat_object@assays$RNA@data

# use mapIds method to obtain Entrez IDs
mapping = mapIds(org.Hs.eg.db, rownames(counts), "ENTREZID", "SYMBOL")
counts <- counts[!is.na(mapping), ]
rownames(counts) <- unname(mapping[!is.na(mapping)])

ccat.v <- CompCCAT(exp = counts, ppiA = net13Jun12.m)
meta$ccat <- ccat.v

write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_MG_SCENT.csv")
