H5PY = c("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/latePRPC/.profile_scores.h5",
"/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/latePRPC/.counts_scores.h5")
MAX_SEQLETS = "50000"
OUTPUT_PREFIX = c("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/latePRPC/profile_scores/",
"/storage/singlecell/zz4/fetal_bash/results/chrombpnet_modisco_motifs/latePRPC/counts_scores/")

sink("cprofile_scores_latePRPC.sh")
cat(paste("chrombpnet modisco_motifs \\
  -i ", H5PY[1], " \\
  -n ",MAX_SEQLETS," \\
  -op ", OUTPUT_PREFIX[1], sep = ""))
cat("\n")
sink()

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_modisco_motifs/counts_scores_latePRPC.sh")
cat(paste("chrombpnet modisco_motifs \\
  -i ", H5PY[2], " \\
  -n ",MAX_SEQLETS," \\
  -op ", OUTPUT_PREFIX[2], sep = ""))
sink()

