dir.create("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/",showWarnings = FALSE)
dir.create("/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_chrom/",showWarnings = FALSE)
rna_input_dir="/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap"
samples = c("AC", "BC", "Cone", "HC", "MG", "RGC", "Rod", "RPC")
adata_rna_file = paste(rna_input_dir,"_",samples,"/adata_umap.h5ad",sep = "")

adata_atac_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/",
                    samples, ".h5ad", sep = "")
output_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/",
                    samples, ".h5ad", sep = "")
df = data.frame(adata_atac_file = adata_atac_file, adata_rna_file= adata_rna_file, samples = samples, output_file = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_chrom/meta.csv",
          row.names = FALSE)