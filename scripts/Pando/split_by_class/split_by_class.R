library(Seurat)
library(Signac)

n_features = 5000

for (x in c("NRPC", "PRPC", "MG")) {
    rna <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_rna/merged_rna.rds")
    meta <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv")
    rownames(meta) <- meta$X
    meta <- meta[meta$majorclass == x, ]
    common_cells <- intersect(colnames(rna), rownames(meta))
    rna <- subset(rna, cells = common_cells)

    atac <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_atac/atac.rds")
    common_cells <- intersect(colnames(rna), colnames(atac))

    rna <- subset(rna, cells = common_cells)
    atac <- subset(atac, cells = common_cells)

    rna <- NormalizeData(rna)
    rna <- ScaleData(rna)
    rna <- FindVariableFeatures(rna, nfeatures = n_features)
    rna <- subset(rna, features = VariableFeatures(rna))

    atac <- RunTFIDF(atac)
    atac <- FindTopFeatures(atac, min.cutoff = "q5")

    atac[["RNA"]] <- rna@assays$RNA

    meta <- meta[common_cells, ]
    atac@meta.data <- cbind(atac@meta.data, meta[colnames(atac), ])

    saveRDS(atac, paste("/storage/chentemp/zz4/adult_dev_compare/results/pando_merged_seurat_object/pando_merged_seurat_object_",
        x, ".rds", sep = ""))
}

for (x in c("RGC", "HC", "Cone", "AC", "BC", "Rod")) {
    rna <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_rna/merged_rna.rds")
    meta <- read.csv(paste("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/",
        x, "_w_NRPC.csv", sep = ""))
    rownames(meta) <- meta$X
    common_cells <- intersect(colnames(rna), rownames(meta))
    rna <- subset(rna, cells = common_cells)

    atac <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/merged_atac/atac.rds")
    common_cells <- intersect(colnames(rna), colnames(atac))

    rna <- subset(rna, cells = common_cells)
    atac <- subset(atac, cells = common_cells)

    rna <- NormalizeData(rna)
    rna <- ScaleData(rna)
    rna <- FindVariableFeatures(rna, nfeatures = n_features)
    rna <- subset(rna, features = VariableFeatures(rna))

    atac <- RunTFIDF(atac)
    atac <- FindTopFeatures(atac, min.cutoff = "q5")

    atac[["RNA"]] <- rna@assays$RNA

    meta <- meta[common_cells, ]
    atac@meta.data <- cbind(atac@meta.data, meta[colnames(atac), ])

    saveRDS(atac, paste("/storage/chentemp/zz4/adult_dev_compare/results/pando_merged_seurat_object/pando_merged_seurat_object_",
        x, ".rds", sep = ""))
}