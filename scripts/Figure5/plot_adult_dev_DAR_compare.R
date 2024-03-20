x_col <- c("BC", "Cone", "Rod", "MG", "NRPC", "PRPC", "RGC", "AC", "HC")
y_col <- c("BC", "Cone", "Rod", "MG", "RGC", "AC", "HC", "DS")
# Specify the number of rows and columns in the matrix
num_rows <- 9
num_cols <- 8

# Create a matrix with all zeros
my_matrix <- matrix(0, nrow = num_rows, ncol = num_cols)

# Set row and column indices
rownames(my_matrix) <- x_col
colnames(my_matrix) <- y_col

for (x in c("BC", "Cone", "Rod", "MG", "NRPC", "PRPC", "RGC", "AC", "HC")){
  for (y in c("BC", "Cone", "Rod", "MG", "RGC", "AC", "HC")) {
    print(x)
    print(y)
    df <- read.csv(paste("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/",
        x, "_Dev_Adult_Compare.bed", sep = ""), sep = "\t", header = FALSE)
    temp <- df[c("V1", "V2", "V3")]
    n1 = nrow(temp[!duplicated(temp), ])
    duplicated_rows <- duplicated(temp) | duplicated(temp, fromLast = TRUE)

    # Get indices of unique rows
    unique_indices <- which(!duplicated_rows)
    df <- df[unique_indices,]
    my_matrix[x,y] = sum(df$V4 == y)
  }
  temp <- read.csv(paste("/storage/chentemp/zz4/adult_dev_compare/results/ArchR/Peak_set/Dev_DAR_peak_",
        x, ".bed", sep = ""), sep = "\t", header = FALSE)
  temp <- temp[!duplicated(temp), ]
  n2 <-nrow(temp)
  my_matrix[x,'DS'] = n2-n1
}

# Load the gplots package
library(gplots)
# Function to create a publication-ready heatmap with Arial font
heatmap_for_publication <- function(data_matrix, col_palette = bluered(100)) {
  # Create a heatmap
  filename = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/plot_adult_dev_DAR_compare.png"
  # Save the heatmap as a publication-quality image if filename is provided
  png(filename, width = 8, height = 6, units = "in", res = 300, bg = "transparent")
  heatmap.2(data_matrix,
            Rowv=FALSE,
            Colv=FALSE,
            cellnote = data_matrix,
            trace="none",
            density.info="none",
            dendrogram = "none",          # Turn off row and column dendrogram
            col = col_palette,       # Set color palette
            scale = "row",           # Scale rows (normalize by row)
            key = FALSE,              # Show color key
            #keysize = 1.5,           # Set the size of the color key
            #density.info = "none",   # Do not show density plot
            margins = c(7, 7),       # Set margins for row and column labels
            #main = "Comparison Between Development DARs and Adult DARs",  # Set the main title
            #xlab = "Adult Major Class",   # Set the x-axis label
            #ylab = "Devlopment Major Class",    # Set the y-axis label
            notecol = 'black',
            cexRow=1.5,
            cexCol = 1.5,
            notecex = 1.5
  )
  dev.off()
    cat("Heatmap saved as", filename, "\n")
  }

heatmap_for_publication(my_matrix)