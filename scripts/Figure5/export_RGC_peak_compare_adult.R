shared <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/RGC_Dev_Adult_Compare.bed", sep = "\t", header = FALSE)
temp <- shared[c("V1", "V2", "V3")]
duplicated_rows <- duplicated(temp) | duplicated(temp, fromLast = TRUE)
unique_indices <- which(!duplicated_rows)
shared <- shared[unique_indices,]
shared <- shared[shared$V4 =="RGC",]

temp <- temp[!duplicated(temp), ]
rownames(temp) <- paste(temp$V1,temp$V2,temp$V3,sep = "_")
dev <- read.csv("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/Dev_DAR_peak_RGC.bed", sep = "\t", header = FALSE)
rownames(dev) <- paste(dev$V1,dev$V2,dev$V3,sep = "_")
dev <- dev[!(rownames(dev) %in% rownames(temp)),]

# Specify the file path where you want to save the BED file
bed_file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/RGC_Dev_Sep.bed"
# Save the dataframe to BED file
write.table(
  dev,
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)

bed_file_path <- "/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/RGC_Shared.bed"
# Save the dataframe to BED file
write.table(
  shared[c("V1", "V2", "V3")],
  bed_file_path,
  sep = "\t",          # Use tab as the column separator
  col.names = FALSE,   # Do not write column names to the file
  quote = FALSE,       # Do not quote the values
  row.names = FALSE    # Do not write row names to the file
)