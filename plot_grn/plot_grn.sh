input_file=(
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/RPC_RPC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/PRPC_PRPC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/NRPC_NRPC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/AC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/Cone_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/Rod_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/HC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/RGC_seurat_object.rds
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/MG_seurat_object.rds
)
output_dir=(
/storage/singlecell/zz4/fetal_bash/figures/MG_Pando/RPC_
/storage/singlecell/zz4/fetal_bash/figures/MG_Pando/PRPC_
/storage/singlecell/zz4/fetal_bash/figures/MG_Pando/NRPC_
/storage/singlecell/zz4/fetal_bash/figures/AC_Pando/
/storage/singlecell/zz4/fetal_bash/figures/BC_Pando/
/storage/singlecell/zz4/fetal_bash/figures/Cone_Pando/
/storage/singlecell/zz4/fetal_bash/figures/Rod_Pando/
/storage/singlecell/zz4/fetal_bash/figures/HC_Pando/
/storage/singlecell/zz4/fetal_bash/figures/RGC_Pando/
/storage/singlecell/zz4/fetal_bash/figures/MG_Pando/
)
for i in "${!input_file[@]}"
do
	slurmtaco.sh -p mhgcp -m 10G --30day -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript plot_grn.R  ${input_file[i]} ${output_dir[i]};
done


