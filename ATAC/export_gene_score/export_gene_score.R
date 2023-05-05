output_results_path="/storage/singlecell/zz4/fetal_snakemake/results/ArchR/gene_score/"
dir.create(output_results_path, showWarnings=FALSE)
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

proj1 <- loadArchRProject("Save-proj1")

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv")
meta$X2 <- stri_replace_last_fixed(meta$X, "_","#")
rownames(meta) <- meta$X2
common_cells <- intersect(rownames(proj1@cellColData),meta$X2)

meta <- meta[common_cells,]
proj2<- subsetArchRProject(proj1,cells = common_cells, outputDirectory = "proj2",force = T)

for (col in colnames(meta)) {
    proj2@cellColData[, col] <- meta[rownames(proj2@cellColData),
        col]
}

gs <- getMatrixFromProject(proj2, useMatrix = "GeneScoreMatrix")
res <- gs@assays@data$GeneScoreMatrix
rownames(res) <- gs@elementMetadata$name

write(colnames(res), file = paste(output_results_path,"gene_score_colnames.txt",sep=""))
write(rownames(res), file = paste(output_results_path,"gene_score_rownames.txt",sep=""))
writeMM(res, file = paste(output_results_path, "gene_score.txt",sep=""))

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-proj2", load = FALSE)