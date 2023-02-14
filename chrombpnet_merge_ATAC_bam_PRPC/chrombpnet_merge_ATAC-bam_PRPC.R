sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_merge_bam/earlyPRPC.sh")
bams <- list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/earlyPRPC",
    recursive = F, pattern = "*.bam$", full.names = T)
cat("samtools merge -f /storage/singlecell/zz4/fetal_bash/results/merged_bam/earlyPRPC_unsorted_ATAC.bam ",
    paste(bams, collapse = " "))
cat("\n")
cat("samtools sort -@4 /storage/singlecell/zz4/fetal_bash/results/merged_bam/earlyPRPC_unsorted_ATAC.bam -o /storage/singlecell/zz4/fetal_bash/results/merged_bam/earlyPRPC_unsorted_ATAC.sort.bam")
cat("\n")
cat("samtools index /storage/singlecell/zz4/fetal_bash/results/merged_bam/earlyPRPC_unsorted_ATAC.sort.bam\n")
sink()

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_merge_bam/latePRPC.sh")
bams <- list.files("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/latePRPC",
    recursive = F, pattern = "*.bam$", full.names = T)
cat("samtools merge -f /storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.bam ",
    paste(bams, collapse = " "))
cat("\n")
cat("samtools sort -@4 /storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.bam -o /storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.sort.bam")
cat("\n")
cat("samtools index /storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.sort.bam\n")
sink()
