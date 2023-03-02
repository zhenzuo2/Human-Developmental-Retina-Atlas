library(ArchR)
library(ComplexHeatmap)
library(circlize)
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(gUtils)
set.seed(0)
seurat_pando_object <- "/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_All_feature_selection_FALSE.rds"
output_dir = "/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC/"
label = "PRPC"
TFs = c("FOXM1", "ESRRG", "NFIB")
pseudotime_file <- "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
set.seed(0)

# Define functions
byapply <- function(x, by, fun, ...) {
    # Create index list
    if (length(by) == 1) {
        nc <- ncol(x)
        split.index <- rep(1:ceiling(nc/by), each = by, length.out = nc)
    } else {
        nc <- length(by)
        split.index <- by
    }
    index.list <- split(seq(from = 1, to = nc), split.index)

    # Pass index list to fun using sapply() and return object
    sapply(index.list, function(i) {
        do.call(fun, list(x[, i], ...))
    })
}

sorting <- function(mat) {
    nrows <- dim(mat)[1]
    ncols <- dim(mat)[2]
    sums <- t(t(mat) * (1:ncols))
    return(mat[order(rowSums(sums)), ])
}

# Read pando result rds object
seurat_object <- readRDS(seurat_pando_object)

# Find and filter modles
seurat_object <- find_modules(seurat_object, p_thresh = 0.1, nvar_thresh = 2,
    min_genes_per_module = 1, rsq_thresh = 0.05)
modules <- NetworkModules(seurat_object)
modules

# Extract tf, target genes and cres
module_info <- data.frame(modules@meta)
cres <- module_info[module_info$tf %in% TFs, "regions"]
n = c()
for (i in 1:length(cres)) {
    n = c(n, stringr::str_count(cres[[i]], ";") + 1)
}
cres <- unlist(stringr::str_split(cres, pattern = ";"))
targets <- rep(module_info[module_info$tf %in% TFs, "target"], n)
tfs <- rep(module_info[module_info$tf %in% TFs, "tf"], n)
chr <- sapply(strsplit(cres, "-"), "[[", 1)
start <- sapply(strsplit(cres, "-"), "[[", 2)
end <- sapply(strsplit(cres, "-"), "[[", 3)

# Read pseudotime
time <- read.csv(pseudotime_file)
rownames(time) <- time$X

# Find nearby peaks
meta <- data.frame(cres = cres, targets = targets, tfs = tfs, chr = chr,
    start = start, end = end)
grange <- makeGRangesFromDataFrame(meta)
dis <- gr.dist(grange, seurat_object@assays$peaks@ranges)
index <- apply(dis, 1, which.min)

# Extract normalized matrix
normalized_matrix_rna <- seurat_object@assays$RNA@data
cells <- intersect(rownames(time), colnames(normalized_matrix_rna))
time <- time[cells, ]
time <- time[order(time$latent_time), ]
cells <- rownames(time)
normalized_matrix_rna <- normalized_matrix_rna[targets, cells]

# Extract pseudotime
latent_time <- time[cells, "latent_time"]

normalized_matrix_atac <- seurat_object@assays$peaks@data
normalized_matrix_atac <- normalized_matrix_atac[index, cells]

normalized_matrix_rna <- byapply(normalized_matrix_rna, 400, rowMeans)
normalized_matrix_atac <- byapply(normalized_matrix_atac, 400, rowMeans)

# weight_rna <- rowSums(normalized_matrix_rna) weight_rna <-
# weight_rna/max(weight_rna) weight_atac <-
# rowSums(normalized_matrix_rna) weight_atac <-
# weight_atac/max(weight_atac)

normalized_matrix_rna <- t(scale(t(normalized_matrix_rna)))
normalized_matrix_atac <- t(scale(t(normalized_matrix_atac)))

# normalized_matrix_rna = apply(t(normalized_matrix_rna), 1,
# function(row) row * weight_rna) normalized_matrix_atac =
# apply(t(normalized_matrix_atac), 1, function(row) row *
# weight_atac)

normalized_matrix_rna <- rbind(rbind(sorting(normalized_matrix_rna[meta$tfs ==
    TFs[1], ]), sorting(normalized_matrix_rna[meta$tfs == TFs[2], ])),
    sorting(normalized_matrix_rna[meta$tfs == TFs[3], ]))
normalized_matrix_atac <- rbind(rbind(sorting(normalized_matrix_atac[meta$tfs ==
    TFs[1], ]), sorting(normalized_matrix_atac[meta$tfs == TFs[2], ])),
    sorting(normalized_matrix_atac[meta$tfs == TFs[3], ]))

svg("/storage/singlecell/zz4/fetal_snakemake/figures/normalized_matrix_rna.svg",
    width = 5, height = 10)
my_cols = viridis::viridis(20, alpha = 1, begin = 0, end = 1, direction = 1, option = "C")

rowlabels = c(rep(TFs[1],sum(meta$tfs == TFs[1])),rep(TFs[2],sum(meta$tfs == TFs[2])),rep(TFs[3],sum(meta$tfs == TFs[3])))
rowlabels = factor(rowlabels, levels = TFs)

Heatmap(as.matrix(normalized_matrix_rna), name = "Z-scores", row_order = 1:(dim(normalized_matrix_rna)[1]),
    #column_split = labels, 
    row_split = rowlabels, 
    #left_annotation = ha2,
    #top_annotation = ha, 
    column_order = 1:(dim(normalized_matrix_rna)[2]), cluster_row_slices = FALSE,
    cluster_column_slices = FALSE, column_names_gp = grid::gpar(fontsize = 0),
    row_names_gp = grid::gpar(fontsize = 0), show_column_names = FALSE,
    col = colorRamp2(seq(-1, 1, length.out = 20), my_cols), row_title = NULL, column_title = NULL,
    use_raster = FALSE)
dev.off()


svg("/storage/singlecell/zz4/fetal_snakemake/figures/normalized_matrix_atac.svg",
    width = 5, height = 10)
my_cols = viridis::viridis(20, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
Heatmap(as.matrix(normalized_matrix_atac), name = "Z-scores", row_order = 1:(dim(normalized_matrix_rna)[1]),
    #column_split = labels, 
    row_split = rowlabels, 
    #left_annotation = ha2,
    #top_annotation = ha, 
    column_order = 1:(dim(normalized_matrix_rna)[2]), cluster_row_slices = FALSE,
    cluster_column_slices = FALSE, column_names_gp = grid::gpar(fontsize = 0),
    row_names_gp = grid::gpar(fontsize = 0), show_column_names = FALSE,
    col = colorRamp2(seq(-1, 1, length.out = 20), my_cols), row_title = NULL, column_title = NULL,
    use_raster = FALSE)
dev.off()