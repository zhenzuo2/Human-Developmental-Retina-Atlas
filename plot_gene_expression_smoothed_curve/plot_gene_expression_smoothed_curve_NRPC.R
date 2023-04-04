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
seurat_pando_object <- "/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/NRPC_seurat_object_All_feature_selection_FALSE.rds"
output_dir = "/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/"
label = "NRPC"
pseudotime_file <- "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/NRPC.obs.csv"
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
TF1 <- c('OTX2', 'NFIB', 'VSX2', 'SMAD2', 'ZEB1', 'RORB', 'LHX4', 'NEUROD1',
       'DACH1', 'RXRG', 'TOX3', 'RREB1', 'NFIA', 'NEUROD4', 'TCF4', 'MYBL1',
       'CRX', 'ZBTB20', 'RAX', 'DACH2', 'FOXP2', 'FOXO3', 'FOXN4', 'NEUROG1',
       'NRL', 'RORA', 'NFIX', 'TCF7L1', 'TCF12', 'RBPJ', 'VSX1', 'RARB',
       'FOXO6', 'ASCL1', 'PRDM8', 'MLLT10', 'INSM1', 'MXI1', 'FOXP4', 'ZBTB7A',
       'NEUROG2', 'TCF7L2', 'STAT3', 'TCF3', 'PPARA', 'HEY2', 'NFIC', 'MEF2C',
       'EGR1')
TF2 <- c('ONECUT1', 'ONECUT2', 'ONECUT3', 'TFAP2B', 'PAX6', 'ESRRG', 'TSHZ2',
       'ZFHX3', 'PROX1', 'ZFPM2', 'ZEB2', 'SOX4', 'TFAP2A', 'EBF1', 'BACH2',
       'SOX5', 'MTA3', 'MEIS1', 'PRDM6', 'NR6A1', 'MYT1L', 'ZNF521', 'FOXN3',
       'ARID5B', 'KLF7', 'PBX1', 'RFX3', 'TSHZ3', 'NHLH2', 'HIVEP3', 'TBX5',
       'ZNF618', 'BCL11A', 'ESRRB', 'ELF2', 'HIVEP2', 'TRERF1', 'NR2F2',
       'HMX1', 'SCRT2', 'DPF1', 'ZNF33B', 'NR4A2', 'ZNF714', 'LHX1', 'ZNF93',
       'ZNF331', 'RFX4', 'TFDP1')
TF3 <- c('ST18', 'HIF1A', 'CCDC88A', 'SOX11', 'PRDM13', 'ZNF423', 'ZNF385D',
       'HEYL', 'ZNF292', 'HES6', 'HIVEP1', 'HES5', 'BARHL2', 'PTF1A', 'TFAP2C',
       'ZBTB18', 'ID2', 'POU2F2', 'TBX3', 'NPAS2', 'PKNOX2', 'FOXO1', 'ZBTB7C',
       'ZNF462', 'HLF', 'SCRT1', 'LCORL', 'ID1', 'NHLH1', 'NFATC3', 'TFDP2',
       'NKX6-1', 'MAF')
TF4 <- c('EBF2', 'PBX3', 'GLI3', 'ATOH7', 'POU4F2', 'HEY1', 'POU6F2', 'ZNF536',
       'LEF1', 'FOXP1', 'HES4', 'SP8', 'LIN28B', 'E2F1', 'NCOA1', 'ISL1',
       'MYB', 'EBF3', 'GLI2', 'CUX2', 'ZBTB16', 'REST', 'SREBF2', 'HMGA2',
       'DLX1', 'TOX', 'ETV5', 'TEAD1', 'SOX2', 'TRPS1', 'TGIF2', 'NPAS3',
       'DLX2', 'ZNF730', 'TFAP2D', 'SMAD7', 'HES1', 'HMGB1', 'PRDM5', 'TOX2',
       'ETV1', 'INSM2', 'E2F5', 'E2F3', 'ID3', 'ZNF519', 'NR2E1', 'PBX4',
       'ZBTB38', 'ZNF732', 'SMAD6', 'ZNF367', 'HMGB2', 'NR3C2', 'ZIC1',
       'MECOM', 'E2F7', 'SMAD3', 'ZNF718', 'NSD2', 'ARNTL2', 'LHX2', 'SOX9',
       'MYBL2', 'ZNF492', 'MITF', 'ZNF257', 'IKZF2', 'HMGXB4', 'E2F2', 'PLAG1',
       'CPEB1', 'TGIF1', 'GLIS3', 'YBX1', 'ZIC5')
TF5 <- c('PRDM1', 'THRB', 'ZHX1', 'SOX6', 'MEIS2', 'ZNF83', 'POU2F1', 'SHOX2',
       'LHX9', 'RUNX1', 'ISL2', 'JUN', 'DDIT3', 'MEF2A', 'HIF3A', 'HMGA1',
       'RXRA', 'RFX2')
TFs <- c(TF1,TF2,TF3,TF4,TF5)
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
time <- time[order(time$velo_s_norm_pseudotime), ]
cells <- rownames(time)
normalized_matrix_rna <- normalized_matrix_rna[unique(c(targets)), cells]

# Extract pseudotime
latent_time <- time[cells, "velo_s_norm_pseudotime"]

normalized_matrix_rna <- byapply(normalized_matrix_rna, 1000, rowMeans)

normalized_matrix_rna <- t(scale(t(normalized_matrix_rna)))

normalized_matrix_rna <- data.frame(normalized_matrix_rna)
normalized_matrix_rna$X <- rownames(normalized_matrix_rna)

df <- normalized_matrix_rna %>%
    pivot_longer(cols = names(normalized_matrix_rna)[1:ncol(normalized_matrix_rna) -
        1], values_to = "newvar")
df$name <- factor(df$name, levels = unique(df$name))
svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/Module1.svg")
ggplot(data = df[df$X %in% TF1, ], aes(x = name, y = newvar, group = X)) +
    geom_line(color = "#1B9E77") + geom_point(color = "#1B9E77") + theme_bw() +
    theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("BC/Rod") +
    ylab("TFs expression") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/Module2.svg")
ggplot(data = df[df$X %in% TF2, ], aes(x = name, y = newvar, group = X)) +
    geom_line(color = "#D95F02") + geom_point(color = "#D95F02") + theme_bw() +
    theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("HC") +
    ylab("TFs expression") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/Module3.svg")
ggplot(data = df[df$X %in% TF3, ], aes(x = name, y = newvar, group = X)) +
    geom_line(color = "#7570B3") + geom_point(color = "#7570B3") + theme_bw() +
    theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("AC") +
    ylab("TFs expression") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/Module4.svg")
ggplot(data = df[df$X %in% TF4, ], aes(x = name, y = newvar, group = X)) +
    geom_line(color = "#E7298A") + geom_point(color = "#E7298A") + theme_bw() +
    theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("Cone") +
    ylab("TFs expression") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

svg("/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/NRPC/Module5.svg")
ggplot(data = df[df$X %in% TF5, ], aes(x = name, y = newvar, group = X)) +
    geom_line(color = "#66A61E") + geom_point(color = "#66A61E") + theme_bw() +
    theme(text = element_text(size = 21)) + xlab("Time") + ggtitle("RGC") +
    ylab("TFs expression") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()
