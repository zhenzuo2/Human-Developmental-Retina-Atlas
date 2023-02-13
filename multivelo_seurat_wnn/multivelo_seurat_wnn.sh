input_file=/storage/singlecell/zz4/fetal_bash/results/Pando_merged/seurat_object.rds
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/MG.csv
)

samples=(
AC
BC
Cone
HC
RGC
Rod
NRPC
PRPC
RPC
MG
)
output_file_path=/storage/singlecell/zz4/fetal_bash/results/multivelo_seurat_wnn/

for i in "${!meta_file[@]}"; do
    slurmtaco.sh -p gpu -m 200G -t 1 -- Rscript multivelo_seurat_wnn.R ${input_file} ${meta_file[i]} ${samples[i]} ${output_file_path};
done
