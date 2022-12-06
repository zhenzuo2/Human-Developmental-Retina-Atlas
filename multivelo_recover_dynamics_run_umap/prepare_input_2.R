df = read.csv( "/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_run_umap/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_run_umap/multivelo_recover_dynamics_run_umap.sh")
for (i in 1:nrow(df)){
  cat("slurmtaco.sh -p gpu -m 50G -t 1 --30day -- python3 multivelo_recover_dynamics_run_umap.py ", df$adata_rna_file[i], " ", df$adata_atac_file[i]," ", df$output_file[i],"\n", sep = "")
}
sink()