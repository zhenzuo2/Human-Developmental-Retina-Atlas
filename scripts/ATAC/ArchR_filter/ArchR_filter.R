output_results_path="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
dir.create("/storage/chentemp/zz4/adult_dev_compare/results/ATAC_filtered_cells",showWarnings = F)
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

cells <- getCellNames(proj1)

write.table(cells, "/storage/chentemp/zz4/adult_dev_compare/results/ATAC_filtered_cells/ATAC_Filtered_cell_before_filterDoublets.csv")

proj2 <- filterDoublets(proj1)

cells <- getCellNames(proj2)

write.table(cells, "/storage/chentemp/zz4/adult_dev_compare/results/ATAC_filtered_cells/ATAC_Filtered_cell.csv")
