files <- c("/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/AC_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/BC_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/Cone_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/HC_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/MG_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/RGC_DAR_peak_hg38",
    "/storage/chentemp/zz4/adult_dev_compare/data/adult_peaks/Rod_DAR_peak_hg38")

labels = c("AC", "BC", "Cone", "HC", "MG", "RGC", "Rod")

for (i in 1:length(files)) {
    temp <- read.csv(files[i], header = FALSE)
    double_list = strsplit(temp$V1, "[:-]")
    df <- as.data.frame(do.call(rbind, double_list))
    bed_file_path = paste("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/Adult_DAR_peak_",
        labels[i], ".bed", sep = "")
    write.table(
        df,
        bed_file_path,
        sep = "\t",          # Use tab as the column separator
        col.names = FALSE,   # Do not write column names to the file
        quote = FALSE,       # Do not quote the values
        row.names = FALSE    # Do not write row names to the file
    )
}
