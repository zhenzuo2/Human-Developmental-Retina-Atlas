input_file=/storage/singlecell/zz4/fetal_bash/results/Pando_merged/seurat_object.rds
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/filtered_cells_MG.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/filtered_cells_RPC.csv
/storage/singlecell/zz4/fetal_bash/results/MultiVelo_filtered_cells/filtered_cells.csv
)
samples=(
AC
BC
Cone
HC
MG
RGC
Rod
RPC
All
)
output_file_path=/storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/

for i in "${!meta_file[@]}"; do
    slurmtaco.sh --g00 -m 150G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript multivelo_seurat_wnn.R ${input_file} ${meta_file[i]} ${samples[i]} ${output_file_path};
done
