sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks/chrombpnet_call_peaks_earlyPRPC.sh")
input_bam = "/storage/singlecell/zz4/fetal_bash/results/merged_bam/earlyPRPC_unsorted_ATAC.sort.cleaned.bam"
output_folder = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/earlyPRPC/"
cat("macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits -p 0.01 -t ",
    input_bam, " --outdir ", output_folder, "\n", sep = "")
sink()

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_call_peaks/chrombpnet_call_peaks_latePRPC.sh")
input_bam = "/storage/singlecell/zz4/fetal_bash/results/merged_bam/latePRPC_unsorted_ATAC.sort.cleaned.bam"
output_folder = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/"
cat("macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits -p 0.01 -t ",
    input_bam, " --outdir ", output_folder, "\n", sep = "")
sink()
