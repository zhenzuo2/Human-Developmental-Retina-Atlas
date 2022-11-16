#!/bin/bash
#SBATCH --job-name=tmp.ymvBU1qhI4
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.ymvBU1qhI4-%j.out
#SBATCH -e out_slurm/tmp.ymvBU1qhI4-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript ArchR_filter.R /storage/singlecell/zz4/fetal_bash/results/ArchR/ /storage/singlecell/zz4/fetal_bash/results/scPred_meta/meta.csv /storage/singlecell/zz4/fetal_bash/results/ArchR/

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0