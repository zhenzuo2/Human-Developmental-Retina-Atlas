library(ArchR)
library(ComplexHeatmap)
library(circlize)
library(Pando)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidyverse)
set.seed(0)
seurat_pando_object <- "/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_seurat_object_All_feature_selection_FALSE.rds"
output_dir = "/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC/"
label = "PRPC"
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
TFs = c("CREB5", "E2F2", "E2F7", "E2F8", "EGR1", "ESRRG", "FOS", "FOXM1", 
"FOXN4", "FOXP2", "GLIS3", "HES6", "HMGB2", "ID2", "MECOM", "MXD3", 
"MYBL1", "NFIA", "NFIB", "NFIX", "NPAS3", "NR4A1", "PROX1", "RARB", 
"SOX5", "SOX6", "TCF4", "TOX", "TOX3", "ZIC1", "ZNF367", "ZNF730")

TF1 <- c("MXD3","HMGB2","FOXM1","E2F8","FOXN4","E2F7","HES6","MYBL1","ZNF367","ZNF730","E2F2")
TF2 <- c("ESRRG","ZIC1","MECOM","NPAS3","TOX")
TF3 <- c("CREB5", "EGR1", "FOS", "FOXP2", "GLIS3", "ID2", 
         "NFIA", "NFIB", "NFIX", "NR4A1", "PROX1", "RARB", "SOX5", 
         "SOX6", "TCF4", "TOX3")

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

# Extract normalized matrix
normalized_matrix_rna <- seurat_object@assays$RNA@data
cells <- intersect(rownames(time), colnames(normalized_matrix_rna))
time <- time[cells, ]
time <- time[order(time$latent_time), ]
cells <- rownames(time)
normalized_matrix_rna <- normalized_matrix_rna[unique(c(targets)), cells]

# Extract pseudotime
latent_time <- time[cells, "latent_time"]

normalized_matrix_rna <- byapply(normalized_matrix_rna, 4000, rowMeans)

normalized_matrix_rna <- t(scale(t(normalized_matrix_rna)))



normalized_matrix_rna <- data.frame(normalized_matrix_rna)
normalized_matrix_rna$X <- rownames(normalized_matrix_rna)

df <- normalized_matrix_rna %>% 
  pivot_longer(cols = names(normalized_matrix_rna)[1:ncol(normalized_matrix_rna)-1],
               values_to = "newvar") 
df$name <- factor(df$name, levels = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"))
svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC/Module1.svg")
ggplot(data=df[df$X %in% TF1,], aes(x=name, y=newvar, group=X)) +
  geom_line(color="red")+
  geom_point(color="red")+ theme_bw() + theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("Module 1")+ylab("TFs expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        )
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC/Module2.svg")
ggplot(data=df[df$X %in% TF2,], aes(x=name, y=newvar, group=X)) +
  geom_line(color="black")+
  geom_point(color="black")+ theme_bw() + theme(text = element_text(size = 21))+ xlab("Time") + ggtitle("Module 2")+ ylab("TFs expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        )
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/PRPC/Module3.svg")
ggplot(data=df[df$X %in% TF3,], aes(x=name, y=newvar, group=X)) +
  geom_line(color="blue")+
  geom_point(color="blue")+ theme_bw() + theme(text = element_text(size = 21))+ xlab("Time") + ggtitle("Module 3")+ ylab("TFs expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        )
dev.off()