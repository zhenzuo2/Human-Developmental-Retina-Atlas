output_results_path = "/storage/singlecell/zz4/fetal_snakemake/results/ArchR/"
dir.create(output_results_path, showWarnings = FALSE)
library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)
setwd("/storage/singlecell/zz4/fetal_snakemake/results/ArchR/")
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 10)

proj2 <- loadArchRProject("Save-proj2")

seRNA <- readRDS("/storage/singlecell/zz4/fetal_snakemake/results/seRNA/seRNA.rds")
proj2 <- addGeneExpressionMatrix(input = proj2, seRNA = seRNA, force = TRUE)

meta <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/cellbygene/adata.obs.csv")
cells <- getCellNames(proj2)
meta$X <- stri_replace_last_fixed(meta$X, "_", "#")
rownames(meta) <- meta$X
common_cells <- intersect(cells, meta$X)
proj2 <- subsetArchRProject(proj2, cells = common_cells, force = TRUE)

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = "proj2")
proj2 <- addReproduciblePeakSet(ArchRProj = proj2, groupBy = "majorclass",
    pathToMacs2 = findMacs2())

proj2 <- addPeakMatrix(proj2)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "PeakMatrix", 
    groupBy = "majorclass",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  labelRows = FALSE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)

proj2 <- addIterativeLSI(ArchRProj = proj2, useMatrix = "TileMatrix", name = "IterativeLSI",
    iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000,
        n.start = 10), varFeatures = 25000, dimsToUse = 1:30)
proj2 <- addPeak2GeneLinks(ArchRProj = proj2, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix")

p <- plotPeak2GeneHeatmap(ArchRProj = proj2, groupBy = "majorclass",k=3,
palATAC = paletteContinuous("horizonExtra"),
  palRNA = paletteContinuous("greenBlue"),)
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp.svg",
    width = 10, height = 8)
p
dev.off()

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-ALL", load = FALSE)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj2,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

df = as.data.frame(table(p2g$idxRNA))
df$Var1 = metadata(p2g)[[2]]$name[df$Var1]
p <- ggplot(df,aes(x=Freq)) + 
  xlab("Number of Enhancers per Gene") + ylab("")+
  geom_bar()+ theme_bw()+theme(text = element_text(size = 30))
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp2.svg",
    width = 10, height = 8)
p
dev.off()

df1 <- as.data.frame(metadata(p2g)[[1]][p2g$idxATAC,])
df2 <- as.data.frame(metadata(p2g)[[2]][p2g$idxRNA,])

colnames(df2) <- c("seqnames_RNA","start_RNA","end_RNA","width_RNA","strand_RNA","name_RNA","idx_RNA")
df <- cbind(df1,df2)
in_or_out <- ((df$start_RNA-df$start)/(df$start_RNA-df$end))<0

#df = data.frame(pmin(abs(df$start_RNA-df$start),abs(df$start_RNA-df$end)))
df = data.frame(abs(df$start_RNA-df$start))

colnames(df)="dis"
df[in_or_out,"dis"] <-0
p <- ggplot(df,aes(x=dis)) + 
  geom_histogram(bins = 100)+
  xlab("Enhancer Gene Pair Distance") + ylab("")+ theme_bw()+theme(text = element_text(size = 30),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp3.svg",
    width = 10, height = 8)
p
dev.off()

atac_meta <- data.frame(metadata(p2g)[[1]])
rownames(atac_meta) <- paste(atac_meta$seqnames,"_",atac_meta$start,"_",atac_meta$end,sep="")
p2g_peaks <- rownames(atac_meta)[p2g$idxATAC]

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
res <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(res) <- c("dars","celltypes")
for (celltypes in c("AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod")){
    dars <- paste(markerList[[celltypes]]$seqnames,"_",markerList[[celltypes]]$start,"_",markerList[[celltypes]]$end,sep="")
    dars <- intersect(dars,p2g_peaks)
    res <- rbind(res,data.frame(dars,celltypes))
}
res <- res[!duplicated(res$dars),]
rownames(res) <- res$dars

data <- data.frame(
  table(res$celltypes)
  )
data <- data[order(data$Freq),]
data$Var1<-factor(data$Var1, levels = data$Var1)
# Barplot
p <- ggplot(data, aes(x=Var1, y=Freq, fill = Var1)) + 
  geom_bar(stat = "identity")+ labs(title = "Cell Type Specific Enhancer Gene Pairs")+
  xlab("") + ylab("Enhancer Gene Pairs")+ theme_bw()+theme(text = element_text(size = 25),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position = "none")
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp4.svg",
    width = 10, height = 10)
p
dev.off()

df1 <- as.data.frame(metadata(p2g)[[1]][p2g$idxATAC,])
df2 <- as.data.frame(metadata(p2g)[[2]][p2g$idxRNA,])

colnames(df2) <- c("seqnames_RNA","start_RNA","end_RNA","width_RNA","strand_RNA","name_RNA","idx_RNA")
df <- cbind(df1,df2)
rownames(df) <- paste(df2$name_RNA,df1$seqnames, df1$start,df1$end,sep = "_")
df$peak <- paste(df1$seqnames, df1$start,df1$end,sep = "_")
in_or_out <- ((df$start_RNA-df$start)/(df$start_RNA-df$end))<0

temp = data.frame(pmin(abs(df$start_RNA-df$start),abs(df$start_RNA-df$end)))
colnames(temp)="dis"
temp[in_or_out,"dis"] <-0

df$dis <- temp$dis
df <- df[df$peak %in% res$dars,]
df$celltypes <- res[df$peak,"celltypes"]
df <- df[df$celltypes %in% c("AC", "BC", "Rod", "Cone", "MG", "HC", "PRPC", "RGC"),]
df$celltypes <- factor(df$celltypes,     # Reorder factor levels
                         c("AC", "BC", "Rod", "Cone", "MG", "HC", "PRPC", "RGC"))

p <- ggplot(df, aes(x=celltypes, y=dis,fill =celltypes )) + 
  geom_boxplot()+ 
  labs(title = "Cell Type Specific Enhancer Gene Pair Distance")+
  xlab("") + ylab("Enhancer Gene Distance")+ theme_bw()+theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position = "none")
p
svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp5.svg",
    width = 10, height = 10)
p
dev.off()


# Load the GenomicRanges package
library(GenomicRanges)
# Function to read genomic ranges
read_genomic_ranges <- function(range_strings) {
  # Initialize an empty GRanges object
  ranges <- GRanges(seqnames = character(0),
                    ranges = IRanges(start = integer(0), end = integer(0)))
  
  # Loop through each range string and add to the GRanges object
  for (range_str in range_strings) {
    range_parts <- unlist(strsplit(range_str, "[:-]"))
    seqname <- range_parts[1]
    start <- as.integer(range_parts[2])
    end <- as.integer(range_parts[3])
    
    new_range <- GRanges(seqnames = seqname,
                         ranges = IRanges(start = start, end = end))
    
    ranges <- c(ranges, new_range)
  }
  
  return(ranges)
}

# Example range strings
range_strings <- c("chr14:56888304-56888662",
                   "chr14:56964627-56965216",
                   "chr14:56902629-56903105",
                   "chr14:56899344-56901557",
                   "chr14:56820328-56822794",
                   "chr14:56807853-56808469",
                   "chr14:56805457-56806309",
                   "chr14:56805457-56807186",
                   "chr14:56730099-56731208",
                   "chr14:56742261-56743173",
                   "chr14:56754441-56756016",
                   "chr14:56735566-56735995",
                   "chr14:56753241-56754311")

# Read the genomic ranges
genomic_ranges <- read_genomic_ranges(range_strings)

# Print the resulting GRanges object
print(genomic_ranges)
markerGenes  <- c("OTX2")
p <- plotBrowserTrack(
    ArchRProj = proj2, 
    groupBy = "majorclass", 
    useGroups = c("BC", "Rod", "Cone","NRPC","AC","MG", "HC", "PRPC", "RGC"),
    geneSymbol = markerGenes, 
    loops = getPeak2GeneLinks(proj2),
    upstream = 150000,
  downstream = 200000,
  features = genomic_ranges
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj2, 
    addDOC = FALSE, width = 8, height = 5)

res <- data.frame(proj2@peakSet)
res$peak_id <- paste(res$seqnames,"_",res$start,"_",res$end,sep="")
atac_meta <- data.frame(metadata(p2g)[[1]])[p2g$idxATAC,]
rownames_atac_meta <- paste(atac_meta$seqnames,"_",atac_meta$start,"_",atac_meta$end,sep="")
res_ <- res[res$peak_id %in% rownames_atac_meta,]

df <- data.frame(values = c(res$distToTSS,res_$distToTSS),id = c(rep("all",nrow(res)),rep("p2g",nrow(res_))))
library(ggplot2)
# Basic density
library(plyr)
df$value2 = 0
mu <- ddply(df, "id", summarise, grp.mean=mean(values))
df[df$values<1000,"values2"] = "0-1"
df[(1000<=df$values)&(df$values<3000),"values2"] = "1-3"
df[(3000<=df$values)&(df$values<5000),"values2"] = "3-5"
df[(5000<=df$values)&(df$values<10000),"values2"] = "5-10"
df[(10000<=df$values)&(df$values<100000),"values2"] = "10-100"
df[100000<=df$values,"values2"] = ">100"
df$value2 <- factor(df$value2,levels = c("0-1" , "1-3","3-5","5-10","10-100",">100"))
# library
library(ggplot2)

# create a dataset
specie <- c(rep("All Peaks" , 6) , rep("Peak-Gene pairs" , 6))
condition <- rep(c(">100","0-1","1-3","3-5","5-10","10-100") , 2)
value <- c(table(df[df$id=='all','values2'])/sum(df$id=='all'),table(df[df$id=='p2g','values2'])/sum(df$id=='p2g'))
data <- data.frame(specie,condition,value)
data$condition <- factor(data$condition,levels = c("0-1" , "1-3","3-5","5-10","10-100",">100"))
# library
# Stacked + percent
p <- ggplot(data, aes(fill=specie, y=value, x=condition)) + theme_bw()+
    geom_bar(stat="identity",position=position_dodge())+xlab("Distance between peak and transcription start site (TSS) in kb")+ylab("Proportion")+ guides(fill=guide_legend(title=""))+theme(text = element_text(size = 20))

svg("/storage/singlecell/zz4/fetal_snakemake/temp/temp6.svg",
    width = 15, height = 10)
p
dev.off()



### Save DARpeak to bed
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
res <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(res) <- c("dars","celltypes","seqnames","start","end")
for (celltypes in c("AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod")){
    dars <- paste(markerList[[celltypes]]$seqnames,"_",markerList[[celltypes]]$start,"_",markerList[[celltypes]]$end,sep="")

    res <- rbind(res,data.frame(dars,celltypes,markerList[[celltypes]]$seqnames,markerList[[celltypes]]$start,markerList[[celltypes]]$end))
}
res <- res[!duplicated(res$dars),]
colnames(res) <- c("dars","celltypes","seqnames","start","end")
rownames(res) <- res$dars

res <- res[res$seqnames != "chrX",]
res$seqnames <- gsub("chr", "", res$seqnames)
# Save DataFrame to BED file using write.table()
write.table(res[,c("seqnames","start","end")], file = "/storage/singlecell/zz4/fetal_snakemake/data/Histone/DAR.bed", sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE)

for (x in unique(res$celltypes)){
  write.table(res[res$celltypes==x,c("seqnames","start","end")], file = paste("/storage/singlecell/zz4/fetal_snakemake/data/Histone/",x,"_DAR.bed",sep=""), sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE)
}