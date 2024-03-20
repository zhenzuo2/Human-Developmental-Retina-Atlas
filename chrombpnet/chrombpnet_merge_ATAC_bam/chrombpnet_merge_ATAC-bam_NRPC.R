merge_clean_bam <- function(input_bam_folder, name, output_sh, output_folder) {
    merged_bam = paste(output_folder, name, "_unsorted_ATAC.bam", sep = "")
    merged_sorted_bam = paste(output_folder, name, "_unsorted_ATAC.sort.bam",
        sep = "")
    merged_sorted_cleaned_bam = paste(output_folder, name, "_unsorted_ATAC.sort.cleaned.bam",
        sep = "")
    sink(output_sh)
    bams <- list.files(paste(input_bam_folder, name, "/", sep = ""), recursive = F,
        pattern = "*.bam$", full.names = T)

    if (!file.exists(merged_bam)) {
        cat("samtools merge -f ", merged_bam, " ", paste(bams, collapse = " "),
            sep = "")
        cat("\n")
    }
    if (!file.exists(merged_sorted_bam)) {
        cat("samtools sort -m 10G -@ 8 ", merged_bam, " -o ", merged_sorted_bam,
            "\n", sep = "")
        cat()
        cat("samtools index ", merged_sorted_bam, "\n", sep = "")
    }
    if (!file.exists(merged_sorted_cleaned_bam)) {
        cat("samtools view -o ", merged_sorted_cleaned_bam, " ", merged_sorted_bam,
            " `seq 1 22 | sed 's/^/chr/'`\n", sep = "")
        cat()
        cat("samtools index ", merged_sorted_cleaned_bam, "\n", sep = "")
    }
    sink()
}

input_bam_folder = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/"
names = c("AC_NRPC", "BC_Rod_NRPC", "Cone_NRPC", "HC_NRPC", "RGC_NRPC")
for (name in names) {
    output_sh = paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_merge_bam/",
    name, ".sh", sep = "")
    output_folder = "/storage/singlecell/zz4/fetal_bash/results/merged_bam/"
    merge_clean_bam(input_bam_folder, name, output_sh, output_folder)
}