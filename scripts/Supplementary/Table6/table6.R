library(Seurat)
library(Pando)
library(tidyr)
library(xlsx)
file = "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table6/Supplementary Table 6.xlsx"
cell_type = "BC"
seurat_object <- readRDS(paste("/storage/chentemp/zz4/adult_dev_compare/results/pando_results/pando_",
    cell_type, ".rds", sep = ""))
seurat_object <- find_modules(seurat_object)
modules <- NetworkModules(seurat_object)
write.xlsx(modules@meta, file, sheetName = paste(cell_type, " GRN"), col.names = TRUE,
    row.names = TRUE, append = F)
for (cell_type in c("AC", "HC", "RGC", "Rod", "Cone", "MG", "PRPC", "NRPC")) {
    seurat_object <- readRDS(paste("/storage/chentemp/zz4/adult_dev_compare/results/pando_results/pando_",
        cell_type, ".rds", sep = ""))
    seurat_object <- find_modules(seurat_object)
    modules <- NetworkModules(seurat_object)
    write.xlsx(as.data.frame(modules@meta), file, sheetName = paste(cell_type, " GRN"),
        col.names = TRUE, row.names = FALSE, append = T)
}