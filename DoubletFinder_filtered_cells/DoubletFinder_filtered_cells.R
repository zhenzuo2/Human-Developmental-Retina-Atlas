meta <- c()
files <- list.files("/storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object",
    full.names = T)
for (f in files) {
    print(f)
    seurat_object <- readRDS(f)
    meta <- c(meta, rownames(seurat_object@meta.data))
}
write.table(meta, file = "/storage/singlecell/zz4/fetal_bash/results/DoubletFinder_filtered_cells/DoubletFinder_filtered_cells.csv",
    row.names = F, col.names = F)
