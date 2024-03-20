library(ArchR)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

input_path = "/storage/singlecell/zz4/fetal_bash/results/ArchR/"
set.seed(1)
setwd(input_path)

PRPC <- loadArchRProject("PRPC")
PRPC$Days <- as.character(PRPC$Days)
PRPC$Days_ <- mapvalues(PRPC$Days, from = c(70, 79, 87, 91, 100, 103, 116,
    136, 137, 141, 142, 162, 165), to = c("FW10", "FW10", "FW13", "FW13",
    "FW13", "FW16", "FW16", "FW19", "FW19", "FW19", "FW19", "FW23", "FW23"),
    warn_missing = TRUE)

save_meta <- function(dat, path) {
    for (sample in unique(dat$sampleid)) {
        print(sample)
        cellid <- rownames(dat[dat$sampleid == sample, ])
        cellid <- str_split_i(cellid, "#", 2)
        write.table(paste("CB:Z:", cellid, sep = ""), file = paste(path,
            sample, ".csv", sep = ""), quote = FALSE, row.names = FALSE,
            col.names = FALSE)
    }
}

save_meta(PRPC@cellColData[PRPC$Days_ %in% c("FW10", "FW13"), ], "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/earlyPRPC/")

save_meta(PRPC@cellColData[PRPC$Days_ %in% c("FW16", "FW19", "FW23"), ],
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/latePRPC/")
