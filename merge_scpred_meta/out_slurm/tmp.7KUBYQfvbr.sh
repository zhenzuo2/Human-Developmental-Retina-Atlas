#!/bin/bash
#SBATCH --job-name=tmp.7KUBYQfvbr
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.7KUBYQfvbr-%j.out
#SBATCH -e out_slurm/tmp.7KUBYQfvbr-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript merge_scpred_meta.R /storage/singlecell/zz4/fetal_bash/results/scPred/ /storage/singlecell/zz4/fetal_bash/results/scPred_meta/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
