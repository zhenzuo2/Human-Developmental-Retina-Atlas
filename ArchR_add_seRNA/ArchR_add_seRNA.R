input_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"
input_meta="/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/merged_raw_filtered_umap_10000_major_sub_class.obs.csv"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"

library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(input_path)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

projretina3 <- loadArchRProject("Save-projretina3")

seRNA <- readRDS("/storage/singlecell/zz4/fetal_bash/results/seRNA/seRNA.rds")
projretina4 <- addGeneExpressionMatrix(input = projretina3, seRNA = seRNA,
                                       force = TRUE)

projretina3
getAvailableMatrices(projretina4)

saveArchRProject(ArchRProj = projretina4, outputDirectory = "Save-projretina4",
                 load = FALSE)
