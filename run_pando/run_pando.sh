input_rna_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
input_atac_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss_filtered.rds"
cell_type=(
BC
AC
Cone
Rod
HC
RGC
RPC
PRPC
NRPC
MG
)
n_features=10000
n_cells=40000
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/RPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/NRPC.csv
/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/MG.csv
)

output_dir=(
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/
/storage/singlecell/zz4/fetal_bash/results/RPC_Pando/
/storage/singlecell/zz4/fetal_bash/results/PRPC_Pando/
/storage/singlecell/zz4/fetal_bash/results/NRPC_Pando/
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "FALSE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "All";
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "FALSE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "Macula";
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "FALSE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "Peripheral";
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "TRUE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "All";
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "TRUE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "Macula";
	slurmtaco.sh --g01 -m 50G --30day -t 8 -- Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} "TRUE" ${n_features} ${n_cells} ${output_dir[i]} "TRUE" "Peripheral";
done