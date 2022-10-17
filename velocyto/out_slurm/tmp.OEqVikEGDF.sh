#!/bin/bash
#SBATCH --job-name=tmp.OEqVikEGDF
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=7-00:00:00
#SBATCH --nodelist=mhgcp-g00
#SBATCH -o out_slurm/tmp.OEqVikEGDF-%j.out
#SBATCH -e out_slurm/tmp.OEqVikEGDF-%j.err

start=$(date +%s)
echo "starting at $(date) on $(hostname)"

# Print the SLURM job ID.
echo "SLURM_JOBID=$SLURM_JOBID"

# Run the application
velocyto run -o /storage/singlecell/zz4/fetal_bash/results/RNA_velocity/Multiome_16w4d_NR -@ 12 --samtools-memory 10240 --bcfile /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_16w4d_NR/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -e Multiome_16w4d_NR /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_16w4d_NR/outs/gex_possorted_bam.bam /storage/singlecell/zz4/Reference/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf

end=$(date +%s)
echo "ended at $(date) on $(hostname). Time elapsed: $(date -u -d @$((end-start)) +'%H:%M:%S')"
exit 0
