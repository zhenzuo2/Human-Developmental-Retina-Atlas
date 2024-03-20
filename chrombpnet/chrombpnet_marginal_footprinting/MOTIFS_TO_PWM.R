TFs_pando <- read.csv("/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/PRPC_modules_meta_All_feature_selection_FALSE.csv")
TFs <- unique(TFs_pando$tf)

early_motifs <- c(list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/earlyPRPC/counts_scores/_reports/",
    pattern = "*.png"), list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/earlyPRPC/profile_scores/_reports/",
    pattern = "*.png"))
early_motifs <- unique(sub("_.*", "", early_motifs))

late_motifs <- c(list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/latePRPC/counts_scores/_reports/",
    pattern = "*.png"), list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/latePRPC/profile_scores/_reports/",
    pattern = "*.png"))
late_motifs <- unique(sub("_.*", "", late_motifs))

early_TFs <- intersect(early_motifs, TFs)

late_TFs <- intersect(late_motifs, TFs)

setdiff(early_TFs, late_TFs)

setdiff(late_TFs, early_TFs)

genes = c(setdiff(early_TFs, late_TFs), setdiff(late_TFs, early_TFs))

genes = c("POU3F3","NFIX_MA1528.1","NFIX_MA0671.1","NFIA","NFIB","E2F8","E2F2_MA0864.2","E2F2_MA0864.1")
sink("/storage/singlecell/zz4/fetal_bash/data/pmf/PRPC.tsv")
for (gene in genes) {
    pfm <- read.table(paste("/storage/singlecell/zz4/fetal_bash/data/pmf/",
        gene, ".pfm", sep = ""), skip = 1)
    pfm <- as.matrix(pfm)
    # Convert the matrix to a named list
    pfm_list <- lapply(1:ncol(pfm), function(i) pfm[, i])
    # Convert the PFM to a motif string
    motif <- c("A", "C", "G", "T")[unlist(lapply(pfm_list, which.max),
        recursive = TRUE)]
    cat(gene)
    cat("\t")
    cat(paste(motif, collapse = ""))
    cat("\n")
}
sink()

genes = c("OTX2_MA0712.1")
sink("/storage/singlecell/zz4/fetal_bash/data/pmf/NRPC.tsv")
for (gene in genes) {
    pfm <- read.table(paste("/storage/singlecell/zz4/fetal_bash/data/pmf/",
        gene, ".pfm", sep = ""), skip = 1)
    pfm <- as.matrix(pfm)
    # Convert the matrix to a named list
    pfm_list <- lapply(1:ncol(pfm), function(i) pfm[, i])
    # Convert the PFM to a motif string
    motif <- c("A", "C", "G", "T")[unlist(lapply(pfm_list, which.max),
        recursive = TRUE)]
    cat(gene)
    cat("\t")
    cat(paste(motif, collapse = ""))
    cat("\n")
}
sink()
