#!/bin/bash
#SBATCH --job-name=tmp.mwToja6hlq
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.mwToja6hlq-%j.out
#SBATCH -e out_slurm/tmp.mwToja6hlq-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript scPred.R /storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/Multi_Fetal_14w5d_FR.rds /storage/singlecell/zz4/fetal_bash/results/scPred/ Multi_Fetal_14w5d_FR /storage/chen/data_share_folder/jinli/scpred/scPred_trainmodel_RNA_svmRadialWeights_scpred.rds /storage/singlecell/zz4/fetal_bash/data/Retina_fetal_sample_meta.csv

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
