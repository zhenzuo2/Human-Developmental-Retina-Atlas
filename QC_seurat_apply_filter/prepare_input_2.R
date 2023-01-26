df = read.csv("/storage/singlecell/zz4/fetal_bash/scripts/QC_seurat_apply_filter/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/QC_seurat_apply_filter/QC_seurat_apply_filter.sh")
cat("input_path=/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/\n")
cat("output_results_path=/storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object/\n")
cat("output_figures_path=/storage/singlecell/zz4/fetal_bash/figures/after_qc_seurat_object/\n")
for (i in 1:nrow(df)) {
  cat(paste("slurmtaco.sh -p short -m 10G -t 1 -- Rscript QC_seurat_apply_filter.R",
            " ", df$Samples_ID[i], " ", df$Samples[i],
            " $output_results_path", " $output_figures_path;", sep = ""))
  cat("\n")
}
sink()