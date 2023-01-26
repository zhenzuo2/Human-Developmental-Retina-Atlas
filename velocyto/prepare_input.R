
# Please see more details at
# https://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples
# Passing parameters from the command line to this script
library(optparse)

option_list = list(make_option(c("-i", "--input"), type = "character",
                               default = "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal",
                               help = "Input folder path name for all samples, for example, /storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"),
                   make_option(c("-o", "--output"), type = "character", default = "/storage/singlecell/zz4/fetal_bash/results/",
                               help = "Output folder path name"), make_option(c("-g", "--gtf"),
                                                                              type = "character", default = "/storage/singlecell/zz4/Reference/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf",
                                                                              help = "genome annotation file, please use the same one used in cellranger for 10X Chromium samples. "))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

samples <- list.dirs(opt$input, full.names = F, recursive = F)

if (substr(opt$input, nchar(opt$input), nchar(opt$input)) != "/") {
  opt$input = paste(opt$input, "/", sep = "")
}
if (substr(opt$output, nchar(opt$output), nchar(opt$output)) != "/") {
  opt$output = paste(opt$output, "/", sep = "")
  dir.create(file.path(opt$output, "RNA_velocity"), showWarnings = FALSE)
}

sink("/storage/singlecell/zz4/fetal_bash/scripts/velocyto/velocyto.sh")
cat("\n")
cat("export LC_ALL=en_US.utf-8")
cat("\n")
cat("export LANG=en_US.utf-8")
cat("\n")

dir.create(paste(opt$output, "RNA_velocity/", sep = ""), showWarnings = FALSE)

for (sam in samples) {
  dir.create(paste(opt$output, "RNA_velocity/", sam, sep = ""), showWarnings = FALSE)
  if (file.exists(paste(opt$input, sam, "/outs/gex_possorted_bam.bam",sep = ""))){
    cat(paste("slurmtaco.sh -p gpu -m 100G -t 1 -- velocyto run -o ", opt$output, "RNA_velocity/", sam, " -@ 12 --samtools-memory 10240 --bcfile ",
              opt$input, sam, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -e ",
              sam, " ", opt$input, sam, "/outs/gex_possorted_bam.bam ", opt$gtf,
              sep = ""))
    cat("\n")
  }
  if (file.exists(paste(opt$input, sam, "/outs/possorted_genome_bam.bam",sep = ""))){
    cat(paste("slurmtaco.sh -p gpu -m 100G -t 1 -- velocyto run -o ", opt$output, "RNA_velocity/", sam, " -@ 12 --samtools-memory 10240 --bcfile ",
              opt$input, sam, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -e ",
              sam, " ", opt$input, sam, "/outs/possorted_genome_bam.bam ", opt$gtf,
              sep = ""))
    cat("\n")
  }
}
sink()



