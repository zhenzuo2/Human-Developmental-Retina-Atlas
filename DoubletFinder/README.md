Apply https://github.com/chris-mcginnis-ucsf/DoubletFinder to the dataset.


mamba activate r
rm -rf out_slurm/
Rscript prepare_input_1.R
Rscript prepare_input_2.R
sh DoubletFinder.sh


