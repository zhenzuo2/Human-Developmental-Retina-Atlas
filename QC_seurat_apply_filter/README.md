Perform QC based on cutoffs from QC_seurat.

conda activate r
rm -rf out_slurm/
Rscript prepare_input_1.R
Rscript prepare_input_2.R
sh QC_seurat_apply_filter.sh

