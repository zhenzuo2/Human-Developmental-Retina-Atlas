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
common_cells <- rownames(proj2@cellColData[(proj2@cellColData$majorclass == "Cone")&(proj2@cellColData$Region !="Whole Eye"),])
proj3 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Region", 
    geneSymbol = c("CNGB3"),
    upstream = 500000,
    downstream = 500000,
    #plotSummary = "featureTrack"
)
plotPDF(p, name = "Plot-Tracks-With-Features-CNGB3", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)


p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Region", 
    geneSymbol = c("PDE6H"),
    upstream = 10000,
    downstream = 10000
)
plotPDF(p, name = "Plot-Tracks-With-Features-PDE6H", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)