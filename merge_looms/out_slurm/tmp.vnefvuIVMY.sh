#!/bin/bash
#SBATCH --job-name=tmp.vnefvuIVMY
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH -o out_slurm/tmp.vnefvuIVMY-%j.out
#SBATCH -e out_slurm/tmp.vnefvuIVMY-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
python3 merge_looms.py

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
