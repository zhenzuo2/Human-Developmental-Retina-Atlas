Add manually annotated cell label to object

rm -rf out_slurm/
conda activate r
Rscript prepare_input.R
conda activate python
sh add_adult_annotated_reference_to_object.sh
