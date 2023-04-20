dir.create("/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/",
    showWarnings = FALSE)
adata_atac_file = "/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/adata_atac.h5ad"

samples = samples = c('AC',
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

nn_idx_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/",
    samples, "_nn_idx.txt", sep = "")
nn_dist_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/",
    samples, "_nn_dist.txt", sep = "")
nn_cells_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/",
    samples, "_nn_cells.txt", sep = "")
output_file = paste("/storage/singlecell/zz4/fetal_bash/results/multivelo_knn_smooth_chrom/",
    samples, ".h5ad", sep = "")

df = data.frame(adata_atac_file = adata_atac_file, samples = samples, nn_idx_file = nn_idx_file,
    nn_dist_file = nn_dist_file, nn_cells_file = nn_cells_file, output_file = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/multivelo_knn_smooth_chrom/meta.csv",
    row.names = FALSE)