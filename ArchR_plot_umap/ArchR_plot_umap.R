library(ArchR)
library(parallel)
library(stringi)
set.seed(1)
setwd(input_path)
library(BSgenome.Hsapiens.UCSC.hg38)
library(svglite)
# prepare to import atac-seq data
addArchRGenome("hg38")
addArchRThreads(threads = 20)

projretina3 <- loadArchRProject("Save-projretina3")
umap <- projretina3@embeddings$UMAP$df
umap$scpred_prediction <- projretina3$scpred_prediction
colnames(umap) <- c("x", "y", "scpred_prediction")
write.csv(umap, "/storage/singlecell/zz4/fetal_bash/results/ArchR/umap.csv")

p <- ggplot(umap, aes(x = x, y = -y, color = scpred_prediction)) + geom_point(size = 0.01) +
    theme_void() + scale_color_manual(values = c(RPC = "#1F77B4", RGC = "#FF7F0E",
    Cone = "#2CA02C", HC = "#D62728", AC = "#9467BD", Rod = "#8C564B",
    BC = "#E377C2", MG = "#7F7F7F")) + theme(plot.background = element_rect(fill = "transparent")) +
    theme(legend.position = "none")
svglite("/storage/singlecell/zz4/fetal_bash/results/ArchR/umap.svg")
p
dev.off()
p
ggsave("/storage/singlecell/zz4/fetal_bash/results/ArchR/umap.png", p, dpi = 600)