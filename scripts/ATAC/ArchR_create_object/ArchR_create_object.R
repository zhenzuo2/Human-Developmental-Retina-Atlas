#https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
#https://www.archrproject.com/

input_meta_file="/storage/chentemp/zz4/adult_dev_compare/scripts/ATAC/ArchR_create_object/meta.csv"
output_results_path="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"

# import packages
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
dir.create(output_results_path, showWarnings = FALSE)
setwd(output_results_path)

# prepare to import atac-seq data
df <- read.csv(input_meta_file)
inputFiles <- df$inputFiles
Samples <- df$Samples
addArchRGenome("hg38")
addArchRThreads(threads = 10)

# remove all unwanted chr
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = Samples,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM", "chrY","KI270728.1", "KI270727.1", "GL000009.2", "GL000194.1", 
                 "GL000205.2", "GL000195.1", "GL000219.1", "KI270734.1", "GL000213.1", 
                 "GL000218.1", "KI270731.1", "KI270721.1", "KI270726.1", "KI270711.1", 
                 "KI270713.1"),subThreading = F, force = T
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# create project
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj1",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

paste0("Memory Size = ", round(object.size(proj1) / 10^6, 3), " MB")

getAvailableMatrices(proj1)

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-proj1", load = FALSE)