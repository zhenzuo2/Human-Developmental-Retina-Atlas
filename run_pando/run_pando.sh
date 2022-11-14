input_rna_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
input_atac_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss_filtered.rds"
cell_type=(
RPC
PRPC
NRPC
AC
BC
Cone
Rod
HC
RGC
MG
)
n_features=5000
meta_file=(
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/PRPC_annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/NRPC_annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap_obs.csv
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_obs.csv
)
output_dir=(
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/RPC_
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/PRPC_
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/NRPC_
/storage/singlecell/zz4/fetal_bash/results/AC_Pando/
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/
/storage/singlecell/zz4/fetal_bash/results/Cone_Pando/
/storage/singlecell/zz4/fetal_bash/results/Rod_Pando/
/storage/singlecell/zz4/fetal_bash/results/HC_Pando/
/storage/singlecell/zz4/fetal_bash/results/RGC_Pando/
/storage/singlecell/zz4/fetal_bash/results/MG_Pando/
)

for i in "${!meta_file[@]}"
do
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "All";
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Macula";
	slurmtaco.sh -p gpu -m 50G --30day -t 4 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_pando.R  ${input_rna_file} ${input_atac_file} ${meta_file[i]} ${cell_type[i]} ${n_features} ${output_dir[i]} "Peripheral";
done


