output_results_path="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
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

meta <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv")
rownames(meta) <- stringi::stri_replace_last_fixed(meta$X,"_","#")
common_cells <- intersect(cells,rownames(meta))

proj3 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)
seRNA <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/seRNA/seRNA.rds")
proj4 <- addGeneExpressionMatrix(input = proj3, seRNA = seRNA,
    force = TRUE)

proj4@cellColData$Time <- meta[rownames(proj4@cellColData),"Time"]
proj4@cellColData$Region <- meta[rownames(proj4@cellColData),"Region"]
proj4@cellColData$Days <- meta[rownames(proj4@cellColData),"Days"]
proj4@cellColData$majorclass <- meta[rownames(proj4@cellColData),"majorclass"]
proj4@cellColData$subclass <- meta[rownames(proj4@cellColData),"subclass"]
proj4@cellColData$celltype <- meta[rownames(proj4@cellColData),"celltype"]

saveArchRProject(ArchRProj = proj4, outputDirectory = "Save-proj2",
                 load = FALSE)