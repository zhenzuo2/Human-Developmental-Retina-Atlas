df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/DoubletFinder/meta.csv")
sink("/storage/singlecell/zz4/fetal_bash/scripts/DoubletFinder/DoubletFinder.sh")
cat("input_path=/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/\n")
cat("output_results_path=/storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/\n")
cat("output_figures_path=/storage/singlecell/zz4/fetal_bash/figures/DoubletFinder_UMAP/\n")

for (i in 1:nrow(df)) {
  cat(paste("slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R",
            " ", df$Samples_ID[i], " ", df$Samples[i],
            " $output_results_path", " $output_figures_path;", sep = ""))
  cat("\n")
}
sink()

