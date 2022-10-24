input_RNA_path="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"
input_ATAC_path="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac.rds"
output_file="/storage/singlecell/zz4/fetal_bash/results/Pando_merged/seurat_object.rds"

slurmtaco.sh -p gpu -m 20G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript Pando_merge_object.R "$input_RNA_path" "$input_ATAC_path" "$output_file";