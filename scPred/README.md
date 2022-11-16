Infer cell fate for filtered cells.


conda activate r
rm -rf out_slurm/
Rscript prepare_input_1.R
Rscript prepare_input_2.R
sh scPred.sh