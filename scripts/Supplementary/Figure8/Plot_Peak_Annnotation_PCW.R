output_results_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
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
proj2$Weeks = plyr::mapvalues(proj2$Days, from = c(59, 70, 76, 79, 87,
    91, 100, 103, 116, 137, 141, 142, 162, 165), to = c("PCW8", "PCW10",
    "PCW10", "PCW10", "PCW13", "PCW13", "PCW15", "PCW15", "PCW15", "PCW19",
    "PCW19", "PCW19", "PCW23", "PCW23"))

proj2$Weeks <- factor(proj2$Weeks, levels = c("PCW8", "PCW10", "PCW13",
    "PCW15", "PCW19", "PCW23"))
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Weeks")

pathToMacs2 <- findMacs2()
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "Weeks", pathToMacs2 = pathToMacs2, maxPeaks = 15000000,peaksPerCell = 50000000)
proj2 <- addPeakMatrix(proj2)

saveArchRProject(ArchRProj = proj2, outputDirectory = "PCW", load = FALSE)
