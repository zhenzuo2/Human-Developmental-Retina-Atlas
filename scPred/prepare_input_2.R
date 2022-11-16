df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/scPred/meta.csv")
sink("/storage/singlecell/zz4/fetal_bash/scripts/scPred/scPred.sh")
cat("meta=/storage/singlecell/zz4/fetal_bash/data/Retina_fetal_sample_meta.csv\n")
cat("reference=/storage/chen/data_share_folder/jinli/scpred/scPred_trainmodel_RNA_svmRadialWeights_scpred.rds\n")
cat("output_results_path=/storage/singlecell/zz4/fetal_bash/results/scPred/\n")
cat("output_figures_path=/storage/singlecell/zz4/fetal_bash/figures/scPred/\n")

for (i in 1:nrow(df)) {
  cat(paste("slurmtaco.sh -p short -m 10G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript scPred.R",
            " ", df$Samples_ID[i], " ", df$Samples[i]," $reference", " $meta",
            " $output_results_path", " $output_figures_path", sep = ""))
  cat("\n")
}
sink()

