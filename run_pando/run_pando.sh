input_rna_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
input_atac_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss_filtered.rds"
cell_type=(
AC
BC
Cone
Rod
HC
RGC
RPC
PRPC
NRPC
MG
)
n_features=5000
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_AC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Cone.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Rod.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_HC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_RGC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/RPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/PRPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC.csv
/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/MG.csv
)
output_dir=(
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/
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
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "All";
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Macula";
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Peripheral";
done

i=1
slurmtaco.sh -p short -m 100G -t 8 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "All";
slurmtaco.sh -p short -m 100G -t 8 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Macula";
slurmtaco.sh -p short -m 100G -t 8 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Peripheral";

