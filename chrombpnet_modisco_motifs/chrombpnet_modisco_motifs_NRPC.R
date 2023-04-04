names = c("AC_NRPC", "BC_Rod_NRPC", "Cone_NRPC", "HC_NRPC", "RGC_NRPC")
for (name in names) {
    H5PY = c(paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/",
        name, "/.profile_scores.h5", sep = ""), paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/",
        name, "/.counts_scores.h5", sep = ""))
    MAX_SEQLETS = "50000"
    OUTPUT_PREFIX = c(paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/",
        name, "/profile_scores/", sep = ""), paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/",
        name, "/counts_scores/", sep = ""))

    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_modisco_motifs/profile_scores_",
        name, ".sh", sep = ""))
    cat("mkdir /storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/",
        name, "\n", sep = "")
    cat("mkdir ", OUTPUT_PREFIX[1], "\n", sep = "")
    cat(paste("chrombpnet modisco_motifs \\
  -i ", H5PY[1], " \\
  -n ",
        MAX_SEQLETS, " \\
  -op ", OUTPUT_PREFIX[1], sep = ""))
    cat("\n")
    sink()

    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_modisco_motifs/counts_scores_",
        name, ".sh", sep = ""))
    cat("mkdir /storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/",
        name, "\n", sep = "")
    cat("mkdir ", OUTPUT_PREFIX[2], "\n", sep = "")
    cat(paste("chrombpnet modisco_motifs \\
  -i ", H5PY[2], " \\
  -n ",
        MAX_SEQLETS, " \\
  -op ", OUTPUT_PREFIX[2], sep = ""))
    sink()
}

