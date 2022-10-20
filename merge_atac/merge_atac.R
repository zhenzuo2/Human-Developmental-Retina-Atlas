# Take all args
args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
plan("multicore", workers = 20)
set.seed(0)

mergeGRangesData <- function(bed_list) {
    peaks = data.frame()
    for (i in 1:length(bed_list)) {
        peaks = rbind(peaks, read.table(file = bed_list[i], col.names = c("chr",
            "start", "end")))
    }
    gr_list = makeGRangesFromDataFrame(peaks)
    combined.peaks <- reduce(x = gr_list)
    # Filter out bad peaks based on length
    combined.peaks@seqinfo@seqlengths <- width(combined.peaks)
    combined.peaks <- combined.peaks[combined.peaks@seqinfo@seqlengths <
        10000 & combined.peaks@seqinfo@seqlengths > 20]
    combined.peaks <- combined.peaks[combined.peaks@seqnames %in% c("chr1",
        "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
        "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")]
    return(combined.peaks)
}

quantify_peaks <- function(meta_list, fragment_list, combined.peaks, sams) {
    seurat_object_list <- NULL
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    # seqlevelsStyle(annotation) <- 'UCSC'
    ucsc.levels <- str_replace(string = paste("chr", seqlevels(annotation),
        sep = ""), pattern = "chrMT", replacement = "chrM")
    seqlevels(annotation) <- ucsc.levels
    genome(annotation) <- "hg38"
    for (i in 1:length(meta_list)) {
        md <- read.table(file = meta_list[i], stringsAsFactors = FALSE,
            sep = ",", header = TRUE, row.names = 1)
        md$gex_barcode <- rownames(md)
        md <- md[(md$is_cell == 1) & (md$excluded_reason == 0), ]
        frags <- CreateFragmentObject(path = fragment_list[i], cells = rownames(md))
        counts <- FeatureMatrix(fragments = frags, features = combined.peaks,
            cells = rownames(md))
        assay <- CreateChromatinAssay(counts, fragments = frags, annotation = annotation)
        seurat_object <- CreateSeuratObject(assay, assay = "peaks")
        # compute LSI seurat_object <- FindTopFeatures(seurat_object)
        # seurat_object <- RunTFIDF(seurat_object) seurat_object <-
        # RunSVD(seurat_object)
        seurat_object <- RenameCells(seurat_object, new.names = paste(sams[i],
            colnames(seurat_object), sep = "_"))
        seurat_object_list[[i]] <- seurat_object
    }
    return(seurat_object_list)
}

bed_list <- c("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR_2/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_NR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_FR/outs/atac_peaks.bed",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_NR/outs/atac_peaks.bed")
meta_list <- c("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR_2/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_NR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_FR/outs/per_barcode_metrics.csv",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_NR/outs/per_barcode_metrics.csv")
fragment_list <- c("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR_2/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_11w2d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_13W_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_14w5d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_19W4d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_20W2d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multi_Fetal_23w1d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_10w_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_12w3d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_14w2d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_16w4d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_20w1d_NR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_FR/outs/atac_fragments.tsv.gz",
    "/storage/singlecell/zz4/fetal_bash/data/Retina_fetal//Multiome_23w4d_NR/outs/atac_fragments.tsv.gz")

sams <- c("Multi_Fetal_11w2d_FR", "Multi_Fetal_11w2d_FR_2", "Multi_Fetal_11w2d_NR",
    "Multi_Fetal_13W_FR", "Multi_Fetal_13W_NR", "Multi_Fetal_14w5d_FR",
    "Multi_Fetal_14w5d_NR", "Multi_Fetal_19W4d_FR", "Multi_Fetal_19W4d_NR",
    "Multi_Fetal_20W2d_FR", "Multi_Fetal_20W2d_NR", "Multi_Fetal_23w1d_FR",
    "Multi_Fetal_23w1d_NR", "Multiome_10w_FR", "Multiome_10w_NR", "Multiome_12w3d_FR",
    "Multiome_12w3d_NR", "Multiome_14w2d_FR", "Multiome_14w2d_NR", "Multiome_16w4d_FR",
    "Multiome_16w4d_NR", "Multiome_20w1d_FR", "Multiome_20w1d_NR", "Multiome_23w4d_FR",
    "Multiome_23w4d_NR")

combined.peaks <- mergeGRangesData(bed_list)
retina <- quantify_peaks(meta_list, fragment_list, combined.peaks, sams)

for (i in 1:length(sams)) {
    retina[[i]]$dataset <- sams[i]
}

seurat_object <- merge(retina[[1]], y = c(retina[2:length(sams)]), project = "retina")
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q0")
seurat_object <- RunSVD(seurat_object)

saveRDS(seurat_object, output_file)