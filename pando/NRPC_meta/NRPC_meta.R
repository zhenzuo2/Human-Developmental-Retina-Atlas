meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.csv")

meta <- meta[meta$leiden %in% c("2", "13", "0", "1", "18", "9", "10", "4",
    "17", "7", "14"), ]

meta$subclass = "Unknown"
clusters = c("2", "13")
meta[meta$leiden %in% clusters, "subclass"] = "AC"
clusters = c("0", "1", "18")
meta[meta$leiden %in% clusters, "subclass"] = "HC"
clusters = c("9", "10")
meta[meta$leiden %in% clusters, "subclass"] = "BC"
clusters = c("4")
meta[meta$leiden %in% clusters, "subclass"] = "Rod"
clusters = c("17")
meta[meta$leiden %in% clusters, "subclass"] = "Cone"
clusters = c("7", "14")
meta[meta$leiden %in% clusters, "subclass"] = "RGC"

table(meta$subclass)

write.csv(meta[meta$subclass %in% c("AC", "HC", "RGC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/AC_HC_RGC.csv")
write.csv(meta[meta$subclass %in% c("AC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/AC.csv")
write.csv(meta[meta$subclass %in% c("HC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/HC.csv")
write.csv(meta[meta$subclass %in% c("RGC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/RGC.csv")

write.csv(meta[meta$subclass %in% c("BC", "Cone", "Rod"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/BC_Rod_Cone.csv")
write.csv(meta[meta$subclass %in% c("BC"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/BC.csv")
write.csv(meta[meta$subclass %in% c("Cone"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Cone.csv")
write.csv(meta[meta$subclass %in% c("Rod"), ], "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Rod.csv")

write.csv(meta, "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC/Annotated_NRPC.csv")
