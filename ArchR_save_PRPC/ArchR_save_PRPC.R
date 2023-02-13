library(ArchR)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

input_path = "/storage/singlecell/zz4/fetal_bash/results/ArchR/"
set.seed(1)
setwd(input_path)

projretina4 <- loadArchRProject("Save-projretina4")

cells <- rownames(projretina4@cellColData)
meta <- read.csv("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC.csv")
rownames(meta) <- meta$X
cells = intersect(cells, stringi::stri_replace_last(meta$X, fixed = "_",
                                                    "#"))


PRPC <- subsetArchRProject(projretina4, cells = cells, force = TRUE)

PRPC <- addIterativeLSI(
  ArchRProj = PRPC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,force = TRUE
)

PRPC <- addClusters(
  input = PRPC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,force = TRUE
)

PRPC <- addUMAP(
  ArchRProj = PRPC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",force = TRUE
)

addArchRThreads(threads = 1)
PRPC <- addGroupCoverages(ArchRProj = PRPC, groupBy = "Clusters",
    threads = 1)
addArchRThreads(threads = 20)

pathToMacs2 <- findMacs2()
PRPC <- addReproduciblePeakSet(
    ArchRProj = PRPC, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

PRPC <- addPeakMatrix(PRPC)

saveArchRProject(ArchRProj = PRPC, outputDirectory = "PRPC",
                 load = FALSE)