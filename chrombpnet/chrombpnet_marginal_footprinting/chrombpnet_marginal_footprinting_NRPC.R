labels = c("AC_NRPC", "BC_Rod_NRPC", "Cone_NRPC", "HC_NRPC")
for (name in labels) {
    MODEL_H5 = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/models/chrombpnet.h5", sep = "")
    REGIONS = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/",
        name, "/auxiliary/filtered.nonpeaks.bed", sep = "")
    GENOME = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
    CHR_FOLD_PATH = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_splits_PRPC/fold_0.json"
    OUTPUT_PREFIX = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_marginal_footprinting/",
        name, "/", sep = "")
    MOTIFS_TO_PWM = "/storage/singlecell/zz4/fetal_bash/data/pmf/NRPC.tsv"
    dir.create(OUTPUT_PREFIX, showWarnings = FALSE)
    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_marginal_footprinting/",
        name, ".sh", sep = ""))
    cat("chrombpnet footprints \\
  -m ", MODEL_H5, " \\
  -r ", REGIONS,
        "\\
  -g ", GENOME, " \\
  -fl ", CHR_FOLD_PATH, " \\
  -op ",
        OUTPUT_PREFIX, "\\
  -pwm_f ", MOTIFS_TO_PWM, sep = "")
    sink()
}