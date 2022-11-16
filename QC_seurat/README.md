Following  
https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
to setup the seurat object for QC.

conda activate r
rm -rf out_slurm/
Rscript prepare_input_1.R
Rscript prepare_input_2.R
sh QC_seurat.sh

