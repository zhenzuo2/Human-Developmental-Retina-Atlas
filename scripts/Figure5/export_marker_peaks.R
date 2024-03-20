library(ArchR)
library(parallel)
library(stringi)
library(dplyr)
set.seed(0)

markersPeaks <- readRDS("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/markersPeaks.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
for (x in markersPeaks@colData@rownames) {
    temp = markerList[[x]]
    csv_file_path = paste("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/Dev_DAR_peak_",
        x, ".csv", sep = "")
    write.table(
        temp,
        sep=",",
        csv_file_path
    )
    bed_file_path = paste("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/Dev_DAR_peak_",
        x, ".bed", sep = "")
    write.table(
        temp[c("seqnames","start","end")],
        bed_file_path,
        sep = "\t",          # Use tab as the column separator
        col.names = FALSE,   # Do not write column names to the file
        quote = FALSE,       # Do not quote the values
        row.names = FALSE    # Do not write row names to the file
    )
}