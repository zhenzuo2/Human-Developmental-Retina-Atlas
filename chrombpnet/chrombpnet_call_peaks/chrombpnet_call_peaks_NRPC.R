names = c("AC_NRPC", "BC_Rod_NRPC", "Cone_NRPC", "HC_NRPC", "RGC_NRPC")
for (name in names) {
    sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks/chrombpnet_call_peaks_",
        name, ".sh", sep = ""))
    input_bam = paste("/storage/singlecell/zz4/fetal_bash/results/merged_bam/",
        name, "_unsorted_ATAC.sort.cleaned.bam", sep = "")
    output_folder = paste("/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/",
        name, "/", sep = "")
    dir.create(output_folder, showWarnings = F)
    cat("macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits -p 0.01 -t ",
        input_bam, " --outdir ", output_folder, "\n", sep = "")
    sink()
}