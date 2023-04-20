output_dir="/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap"
samples = c('AC',
'BC',
'Cone',
'HC',
'RGC',
'Rod',
'NRPC',
'PRPC',
"RPC",
"MG",
"BC_Rod_Cone",
"AC_HC_RGC",
"PRPC_MG")
adata_atac_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/",samples,".h5ad",sep = "")
adata_rna_file = "/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/adata_rna.h5ad"
output_file = paste(output_dir,"_",samples,"/adata_umap.h5ad",sep = "")

df = data.frame(adata_atac_file = adata_atac_file, adata_rna_file= adata_rna_file, samples = samples, output_file = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/multivelo_recover_dynamics_run_umap/meta.csv",
          row.names = FALSE)