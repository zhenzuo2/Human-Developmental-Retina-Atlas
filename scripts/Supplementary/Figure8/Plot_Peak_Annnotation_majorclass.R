output_results_path <-"/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(output_results_path)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj3")

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "majorclass",force = TRUE)

pathToMacs2 <- findMacs2()
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "majorclass", pathToMacs2 = pathToMacs2, maxPeaks = 15000000,peaksPerCell = 500000)
proj2 <- addPeakMatrix(proj2)

saveArchRProject(ArchRProj = proj2, outputDirectory = "majorclass",
                 load = FALSE)