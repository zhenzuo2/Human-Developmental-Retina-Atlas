args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_results_path <- args[2]

library(ArchR)
library(parallel)
library(stringi)
library(MASS)
library(Matrix)

set.seed(1)
setwd(input_path)


projretina3 <- loadArchRProject("Save-projretina3")

gs <- getMatrixFromProject(projretina3, useMatrix = "GeneScoreMatrix")

res <- gs@assays@data$GeneScoreMatrix

rownames(res) <- gs@elementMetadata$name

write(colnames(res), file = paste(output_results_path,"gene_score_colnames.txt",sep=""))
write(rownames(res), file = paste(output_results_path,"gene_score_rownames.txt",sep=""))
writeMM(res, file = paste(output_results_path, "gene_score.txt",sep=""))