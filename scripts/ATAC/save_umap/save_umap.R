output_results_path="/storage/chentemp/zz4/adult_dev_compare/results/ArchR/"
library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(output_results_path)
library(BSgenome.Hsapiens.UCSC.hg38)

# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 1)

proj1 <- loadArchRProject("Save-proj3")

res <- proj1@embeddings$UMAPHarmony$df
res$majorclass <- proj1$majorclass
res$subclass <- proj1$subclass
write.csv(res,"/storage/chentemp/zz4/adult_dev_compare/results/ATAC_umap/ATAC_UMAPHarmony_cor.csv")

# OUTPUT peak set for JUN WANG
df <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/peaks/peakSet.csv")
df <- df[,c("seqnames","start","end")]
write.table(df,"/storage/chentemp/zz4/adult_dev_compare/results/peaks/peakSet.bed",row.names = F,col.names = F,sep = "\t",quote = F)

peaks <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g_metadata_peakSet.csv")
df <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g.csv")
res1 <- peaks[unique(df$idxATAC),c("seqnames","start","end")]
write.table(res1,"/storage/chentemp/zz4/adult_dev_compare/results/peaks/p2g_peakSet.bed",row.names = F,col.names = F,sep = "\t",quote = F)

cA <- getCoAccessibility(
    ArchRProj = proj1,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)
df <- data.frame(metadata(cA)[[1]])
res2 <- df[unique(c(cA$queryHits,cA$subjectHits)),c("seqnames","start","end")]
write.table(res2,"/storage/chentemp/zz4/adult_dev_compare/results/peaks/cA_peakSet.bed",row.names = F,col.names = F,sep = "\t",quote = F)

peak_common <- generics::intersect(res1, res2)  # Apply intersect function
write.table(peak_common,"/storage/chentemp/zz4/adult_dev_compare/results/peaks/peak_common.bed",row.names = F,col.names = F,sep = "\t",quote = F)
