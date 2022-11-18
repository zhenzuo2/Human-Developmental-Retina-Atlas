Merge all metafiles.

conda activate r
rm -rf out_slurm/
Rscript prepare_input.R
sh merge_scpred_meta.sh

