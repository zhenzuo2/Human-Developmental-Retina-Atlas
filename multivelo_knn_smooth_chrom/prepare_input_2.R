df = read.csv("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_knn_smooth_chrom/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_knn_smooth_chrom/multivelo_knn_smooth_chrom.sh")
for (i in 1:nrow(df)){
  cat("slurmtaco.sh -p short -m 50G -t 1 -- python3 multivelo_knn_smooth_chrom.py ", df$adata_atac_file[i], " ", df$samples[i]," ", 
      df$nn_idx_file[i], " ", df$nn_dist_file[i], " ", df$nn_cells_file[i], " ", df$output_file[i],"\n", sep = "")
}
sink()