output_results_path="/storage/singlecell/zz4/fetal_snakemake/results/ArchR/"
dir.create("/storage/singlecell/zz4/fetal_snakemake/results/ATAC_filtered_cells",showWarnings = F)
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(output_results_path)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

proj1 <- loadArchRProject("Save-proj1")

proj2 <- filterDoublets(proj1)

cells <- getCellNames(proj2)

write.table(cells, "/storage/singlecell/zz4/fetal_snakemake/results/ATAC_filtered_cells/ATAC_Filtered_cell.csv")
