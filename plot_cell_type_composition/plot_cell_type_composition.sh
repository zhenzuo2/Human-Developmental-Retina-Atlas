input_file_path=/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/merged_h5ad_adult_annotated_obs.csv
output_file_path=/storage/singlecell/zz4/fetal_bash/figures/cell_composition_stacked_bar_plot/

slurmtaco.sh -p short -m 1G -t 1 -- python3 plot_cell_type_composition.py ${input_file_path} ${output_file_path} "Macula";
slurmtaco.sh -p short -m 1G -t 1 -- python3 plot_cell_type_composition.py ${input_file_path} ${output_file_path} "Peripheral";