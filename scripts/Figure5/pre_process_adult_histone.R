H3K27ac <- read.csv("/storage/chentemp/zz4/adult_dev_compare/data/Histone/H3K27ac",sep = "\t",header=FALSE)
chrs <- c("10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
"1", "20", "21", "22", "2", "3", "4", "5", "6", "7", "8", "9","X", "Y")
H3K27ac <- H3K27ac[H3K27ac$V1 %in% chrs,c("V1","V2","V3")]
H3K27ac<- H3K27ac[!duplicated(H3K27ac), ]
H3K27ac$V1 = paste("chr",H3K27ac$V1,sep="")
bed_file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/H3K27ac.bed"
# Save the dataframe to BED file
write.table(
  H3K27ac,
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)

H3K4me2 <- read.csv("/storage/chentemp/zz4/adult_dev_compare/data/Histone/Hu14_17_18_ret_H3K4me2.narrowPeak",sep = "\t",header=FALSE)
H3K4me2 <- H3K4me2[H3K4me2$V1 %in% chrs,c("V1","V2","V3")]
H3K4me2<- H3K4me2[!duplicated(H3K4me2), ]
H3K4me2$V1 = paste("chr",H3K4me2$V1,sep="")
bed_file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/H3K4me2.bed"
# Save the dataframe to BED file
write.table(
  H3K4me2,
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)