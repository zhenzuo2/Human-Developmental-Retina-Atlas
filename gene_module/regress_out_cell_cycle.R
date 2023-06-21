library(Seurat)
library(Matrix)
plan("multisession", workers = 6)
plan()
options(future.globals.maxSize = 500 * 1024^3)

seurat_object <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/merged_rna/merged_rna.rds")
meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_latent_time.csv")
seurat_object <- subset(seurat_object, cells = meta$X)
data(cc.genes)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded
# with Seurat.  We can segregate this list into markers of G2/M phase
# and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",
    nfeatures = 10000)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))

seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes,
    g2m.features = g2m.genes, set.ident = TRUE)

seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
seurat_object <- ScaleData(seurat_object, vars.to.regress = "CC.Difference",
    features = rownames(seurat_object))
res <- seurat_object@assays$RNA@scale.data
saveRDS(seurat_object, "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/merged_rna_PRPC.CC.Difference.rds")

output_results_path = "/storage/singlecell/zz4/fetal_snakemake/results/merged_rna_PRPC_cell_cycle/"
write(colnames(res), file = paste(output_results_path, "cell_cycle_normed_colnames.txt",
    sep = ""))
write(rownames(res), file = paste(output_results_path, "cell_cycle_normed_rownames.txt",
    sep = ""))
write.table(res, file = paste(output_results_path, "cell_cycle_normed.txt",
    sep = ""), row.names = FALSE, col.names = FALSE)
