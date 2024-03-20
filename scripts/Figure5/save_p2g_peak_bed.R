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
addArchRThreads(threads = 10)

proj <- loadArchRProject("Save-proj3")

p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)
peak_meta <- metadata(p2g)[[1]]
gene_meta <- metadata(p2g)[[2]]

peak_ranges <- data.frame(peak_meta@ranges)
peak_chr <- data.frame(peak_meta@seqnames)

p2g$peak_chr <- peak_chr$peak_meta.seqnames[p2g$idxATAC]
p2g$peak_start <- peak_ranges$start[p2g$idxATAC]
p2g$peak_end <- peak_ranges$end[p2g$idxATAC]

gene_ranges <- data.frame(gene_meta@ranges)
gene_chr <- data.frame(gene_meta@seqnames)

p2g$gene_chr <- gene_chr$gene_meta.seqnames[p2g$idxRNA]
p2g$gene_start <- gene_ranges$start[p2g$idxRNA]
p2g$gene_end <- gene_ranges$end[p2g$idxRNA]

write.csv(p2g,"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/p2g.csv")

p2g <- p2g[c("peak_chr", "peak_start", "peak_end")]
p2g <- p2g[!duplicated(p2g), ]
# Specify the file path where you want to save the BED file
bed_file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/p2g_unique.bed"

# Save the dataframe to BED file
write.table(
  p2g,
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)