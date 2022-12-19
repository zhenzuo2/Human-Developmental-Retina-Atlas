# http://bioconductor.org/packages/release/bioc/vignettes/tricycle/inst/doc/tricycle.html

# load packages
library(tricycle)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(viridis)
# read data and subset progenitors from the whole dataset
seurat_object <- readRDS("/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds")
# set default assay to rna so more features will be included for cell
# cycle prediction later
meta <- read.csv("/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/RPC.csv")
seurat_object <- subset(seurat_object, cells = meta$X)
seurat_object <- NormalizeData(seurat_object)

DefaultAssay(seurat_object) <- "RNA"

# convert seurat object to sce so that it can be used as input for
# tricycle
sce <- as.SingleCellExperiment(seurat_object)
sce <- project_cycle_space(sce, exprs_values = "logcounts", species = "human",
                           gname.type = "SYMBOL")
sce <- estimate_cycle_position(sce)
# assign cell cycle stages using Schwabe method
sce <- estimate_Schwabe_stage(sce, gname.type = "SYMBOL", species = "human")
# save the results to be loaded later
saveRDS(sce, "/storage/singlecell/zz4/fetal_bash/results/cell_cycle/tricycle_sce.rds")
write.csv(sce@colData, "/storage/singlecell/zz4/fetal_bash/results/cell_cycle/cell_cycle_meta.csv")