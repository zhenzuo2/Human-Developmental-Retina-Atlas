input_path="/storage/singlecell/zz4/fetal_snakemake/data/Retina_fetal/"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"

slurmtaco.sh -p short -m 20G -t 10 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript creat_object.R "$input_path" "$output_results_path";