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

proj1 <- loadArchRProject("Save-proj3")

getGroupBW(
  ArchRProj = proj1,
  groupBy = "majorclass",
  normMethod = "None", #"ReadsInTSS",
  tileSize = 100,
  maxCells = 10000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)