output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/gene_score/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 1)

proj1 <- loadArchRProject("Save-proj3")

gs <- getMatrixFromProject(proj1, useMatrix = "GeneScoreMatrix")
res <- gs@assays@data$GeneScoreMatrix
rownames(res) <- gs@elementMetadata$name

write(colnames(res), file = paste(output_results_path, "gene_score_colnames.txt",
    sep = ""))
write(rownames(res), file = paste(output_results_path, "gene_score_rownames.txt",
    sep = ""))
writeMM(res, file = paste(output_results_path, "gene_score.txt", sep = ""))
