df = read.csv("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_chrom/meta.csv")

sink("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_chrom/multivelo_recover_dynamics_chrom.sh")
for (i in 1:nrow(df)){
  cat("slurmtaco.sh --g00 -m 50G -t 8 --30day -- python3 multivelo_recover_dynamics_chrom.py ", df$adata_rna_file[i], " ", df$adata_atac_file[i]," ", df$output_file[i],"\n", sep = "")
}
sink()