meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.csv")

write.csv(meta[meta$subclass %in% c("AC", "HC", "RGC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/AC_HC_RGC.csv")
write.csv(meta[meta$subclass %in% c("AC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/AC.csv")
write.csv(meta[meta$subclass %in% c("HC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/HC.csv")
write.csv(meta[meta$subclass %in% c("RGC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/RGC.csv")

write.csv(meta[meta$subclass %in% c("BC", "Cone", "Rod"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/BC_Rod_Cone.csv")
write.csv(meta[meta$subclass %in% c("BC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/BC.csv")
write.csv(meta[meta$subclass %in% c("Cone"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Cone.csv")
write.csv(meta[meta$subclass %in% c("Rod"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Rod.csv")

write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Annotated_NRPC.csv")

