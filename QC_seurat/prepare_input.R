input_path = "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"
df <- data.frame(Samples_ID = list.files(input_path), Samples = paste(input_path,
    list.files(input_path), "/outs/filtered_feature_bc_matrix.h5", sep = ""))
df <- df[df$Samples_ID %in% c("Multi_Fetal_11w2d_FR", "Multi_Fetal_11w2d_FR_2",
    "Multi_Fetal_11w2d_NR", "Multi_Fetal_13W_FR", "Multi_Fetal_13W_NR",
    "Multi_Fetal_14w5d_FR", "Multi_Fetal_14w5d_NR", "Multi_Fetal_19W4d_FR",
    "Multi_Fetal_19W4d_NR", "Multi_Fetal_20W2d_FR", "Multi_Fetal_20W2d_NR",
    "Multi_Fetal_23w1d_FR", "Multi_Fetal_23w1d_NR", "Multiome_10w_FR",
    "Multiome_10w_NR", "Multiome_12w3d_FR", "Multiome_12w3d_NR", "Multiome_14w2d_FR",
    "Multiome_14w2d_NR", "Multiome_16w4d_FR", "Multiome_16w4d_NR", "Multiome_20w1d_FR",
    "Multiome_20w1d_NR", "Multiome_23w4d_FR", "Multiome_23w4d_NR"), ]
rownames(df) <- NULL
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/QC_seurat/meta.csv",
    row.names = F)

sink("/storage/singlecell/zz4/fetal_bash/scripts/QC_seurat/QC_seurat.sh")
cat("min_cells=10\n")
cat("min_features=200\n")
cat("output_figures_path=/storage/singlecell/zz4/fetal_bash/figures/qc_seurat_object/\n")
cat("output_results_path=/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/\n")
for (i in 1:nrow(df)) {
    cat(paste("slurmtaco.sh -p short -m 20G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript QC_seurat.R",
        " ", df$Samples_ID[i], " ", df$Samples[i], " $min_cells", " $min_features",
        " $output_figures_path", " $output_results_path;", sep = ""))
    cat("\n")
}
sink()