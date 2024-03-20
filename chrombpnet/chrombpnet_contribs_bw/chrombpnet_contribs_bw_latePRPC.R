MODEL_H5 = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_trained_models/latePRPC/models/chrombpnet.h5"
REGIONS = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_peaks/latePRPC/NA_peaks.narrowPeak_no_blacklist.bed"
GENOME = "/storage/singlecell/zz4/Reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
CHROM_SIZES = "/storage/singlecell/zz4/Reference/hg38.chrom.sizes"
OUTPUT_PREFIX = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_contribs_bw/latePRPC/"

sink("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_contribs_bw/latePRPC.sh")
cat(paste("chrombpnet contribs_bw \\
  -m ", MODEL_H5, " \\
  -r ", REGIONS,
    " \\
  -g ", GENOME, " \\
  -c ", CHROM_SIZES, " \\
  -op ", OUTPUT_PREFIX,
    " \\
", sep = ""))
sink()

# export CUDA_VISIBLE_DEVICES=1
# slurmtaco.sh --g01 -m 10G --30day -t 1 -- sh latePRPC.sh