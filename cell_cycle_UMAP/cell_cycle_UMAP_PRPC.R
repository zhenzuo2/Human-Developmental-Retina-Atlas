library(Seurat)

plan("multisession", workers = 6)
plan()
options(future.globals.maxSize= 500*1024^3)

seurat_object <- readRDS("/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna_PRPC.rds")

ifnb.list <- SplitObject(ifnb, split.by = "stim")
data(cc.genes)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))

seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
seurat_object <- ScaleData(seurat_object, vars.to.regress = "CC.Difference", features = rownames(seurat_object))

seurat_object <- FindNeighbors(seurat_object)
seurat_object <- FindClusters(seurat_object)
seurat_object <- RunUMAP(seurat_object,dims = 1:30)
saveRDS(seurat_object,"/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna_PRPC.CC.Difference.rds")
write.csv(seurat_object@reductions$umap@cell.embeddings,"/storage/singlecell/zz4/fetal_bash/results/PRPC_UMAP/PRPC_UMAP.csv")

# cell cycle effects strongly mitigated in PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object), nfeatures.print = 10)