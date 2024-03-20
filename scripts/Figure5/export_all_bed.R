library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)
output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)

setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")

proj <- loadArchRProject("Save-proj3")
res <- data.frame(proj@peakSet)
bed_file_path="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/all_peaks.bed"

write.table(
  res[c("seqnames", "start", "end")],
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)

write.csv(res, "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/all_peaks.csv")