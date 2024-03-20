library(ArchR)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

input_path = "/storage/singlecell/zz4/fetal_bash/results/ArchR/"
set.seed(1)
setwd(input_path)

PRPC <- loadArchRProject("PRPC")
PRPC
Time <- read.csv("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv")
Time$X <- stringi::stri_replace_last(Time$X, fixed = "_", "#")
rownames(Time) <- Time$X

common_cells <- intersect(PRPC$cellNames, Time$X)

Time <- Time[common_cells, ]

save_meta <- function(dat, path) {
    dir.create(path, showWarnings = FALSE)
    for (sample in unique(dat$sampleid)) {
        print(sample)
        cellid <- rownames(dat[dat$sampleid == sample, ])
        cellid <- str_split_i(cellid, "#", 2)
        write.table(paste("CB:Z:", cellid, sep = ""), file = paste(path,
            sample, ".csv", sep = ""), quote = FALSE, row.names = FALSE,
            col.names = FALSE)
    }
}

save_meta(Time[Time$latent_time <= quantile(Time$latent_time, 0.25), ],
    "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC_1/")
save_meta(Time[(Time$latent_time >= quantile(Time$latent_time, 0.25)) &
    (Time$latent_time <= quantile(Time$latent_time, 0.5)), ], "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC_2/")
save_meta(Time[(Time$latent_time >= quantile(Time$latent_time, 0.5)) &
    (Time$latent_time <= quantile(Time$latent_time, 0.75)), ], "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC_3/")
save_meta(Time[(Time$latent_time >= quantile(Time$latent_time, 0.75)),
    ], "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/PRPC_4/")