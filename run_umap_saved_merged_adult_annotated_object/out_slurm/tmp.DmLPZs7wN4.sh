#!/bin/bash
#SBATCH --job-name=tmp.DmLPZs7wN4
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH -o out_slurm/tmp.DmLPZs7wN4-%j.out
#SBATCH -e out_slurm/tmp.DmLPZs7wN4-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python3 run_umap_saved_merged_adult_annotated_object.py /storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated.h5ad /storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
