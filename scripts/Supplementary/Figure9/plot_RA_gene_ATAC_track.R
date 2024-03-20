output_results_path = "/storage/chentemp/zz4/adult_dev_compare/results/ArchR"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
set.seed(0)
setwd("/storage/chentemp/zz4/adult_dev_compare/results/ArchR")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj3")
common_cells <- rownames(proj2@cellColData[(proj2@cellColData$majorclass == "PRPC"),])
proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
proj2 <- addBgdPeaks(proj2)
proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "Motif",
  force = TRUE
)

motifPositions <- getPositions(proj2)
motifs <- c(
    "TBX5","VAX2","ALDH1A1","ALDH1A2","ALDH1A3","CYP26A1","CYP26B1","CYP26C1"
  )
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Region")


seFoot <- getFootprints(
  ArchRProj = proj2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Region",
  flank = 200
)
proj2$Region <- factor(proj2$Region, levels = c("Macula","Peripheral","Whole Eye"))
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj2, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 10,
  flank = 200,
  pal = c("#F8766D",  "#00BFC4","#ffdd05")
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj2, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 10,
  flank = 200,
  pal = c("#F8766D",  "#00BFC4","#ffdd05")
)

p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "Region", 
    geneSymbol = motifs, 
    upstream = 20000,
    downstream = 20000,
    loops = getCoAccessibility(proj2),
    sizes = c(12, 1.5, 3, 4),
    ylim = c(0, 1),
    pal = c("#FF0000", "#E64D00", "#CC9900", "#66B2FF", "#3366FF", "#0000FF")
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-proj4_Region_Motif", load = FALSE)