meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv")
meta <- meta[(meta$majorclass %in% c("PRPC", "MG")) & (meta$Time != "Adult"),
    ]
write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/PRPC_MG.csv")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv")
meta <- meta[(meta$majorclass %in% c("NRPC")) & (meta$Time != "Adult"),
    ]
write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/NRPC.csv")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv")
meta <- meta[(meta$majorclass %in% c("PRPC", "MG", "NRPC")) & (meta$Time !=
    "Adult"), ]
write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/PRPC_MG_NRPC.csv")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv")
meta <- meta[(meta$majorclass %in% c("PRPC")) & (meta$Time !=
    "Adult"), ]
write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/PRPC.csv")