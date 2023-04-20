mkdir /storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/
slurmtaco.sh -p short -m 20G -t 10 -- python3 gene_module_analysis_PRPC.py;