save_meta <- function(dat, path) {
  
  for (sample in unique(dat$sampleid)) {
    dir.create(path, showWarnings = F)
    print(sample)
    cellid <- rownames(dat[dat$sampleid == sample, ])
    print(length(cellid))
    cellid <- str_split_i(cellid, "#", 2)
    write.table(paste("CB:Z:", cellid, sep = ""), file = paste(path,
                                                               sample, ".csv", sep = ""), quote = FALSE, row.names = FALSE,
                col.names = FALSE)
  }
}

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC.csv")
df <- df[df$majorclass == "NRPC", ]
rownames(df) <- df$X
rownames(df) <- stringi::stri_replace_last(rownames(df), fixed = "_", "#")
save_meta(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/AC_NRPC/")
table(df$majorclass)
# 3506

df1 <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC.csv")
df1 <- df1[df1$majorclass == "NRPC", ]
table(df1$majorclass)
# 1168
df2 <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Rod.csv")
df2 <- df2[df2$majorclass == "NRPC", ]
table(df2$majorclass)
# 580
df <- rbind(df1, df2)
df <- df %>%
  dplyr::distinct(X, .keep_all = TRUE)
rownames(df) <- df$X
rownames(df) <- stringi::stri_replace_last(rownames(df), fixed = "_", "#")
save_meta(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/BC_Rod_NRPC/")
table(df$majorclass)

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Cone.csv")
df <- df[df$majorclass == "NRPC", ]
rownames(df) <- df$X
rownames(df) <- stringi::stri_replace_last(rownames(df), fixed = "_", "#")
save_meta(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/Cone_NRPC/")
table(df$majorclass)
# 345

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_HC.csv")
df <- df[df$majorclass == "NRPC", ]
rownames(df) <- df$X
rownames(df) <- stringi::stri_replace_last(rownames(df), fixed = "_", "#")
save_meta(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/HC_NRPC/")
table(df$majorclass)
# 1652

df <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_RGC.csv")
df <- df[df$majorclass == "NRPC", ]
rownames(df) <- df$X
rownames(df) <- stringi::stri_replace_last(rownames(df), fixed = "_", "#")
save_meta(df, "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RGC_NRPC/")
table(df$majorclass)
# 1274


