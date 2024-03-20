# Read Adult TF data
library(ggplot2)
library(tidyr)
library(readxl)
Adult_TF <- read_excel("Desktop/Adult_TF.xlsx", sheet = "Major class top20 RSS TFs")
Adult_TF <- as.data.frame(Adult_TF)
extract_before_underscore <- function(input_vector) {
  result <- gsub("^(.*?)_.*$", "\\1", input_vector)
  return(result)
}

for (index in colnames(Adult_TF)[2:21]) {
  Adult_TF[, index] <- extract_before_underscore(unlist(Adult_TF[, index]))
}

# Read Fetal Data
temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "AC GRN")
fetal_AC <- unique(temp$tf)

temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "BC GRN")
fetal_BC <- unique(temp$tf)

temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "HC GRN")
fetal_HC <- unique(temp$tf)

temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "RGC GRN")
fetal_RGC <- unique(temp$tf)

temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "Rod GRN")
fetal_Rod <- unique(temp$tf)

temp <- read_excel("My Drive/Paper Writing Group/Retina Development paper/Supplementary/sup data/Supplementary Table 4.xlsx",
                   sheet = "Cone GRN")
fetal_Cone <- unique(temp$tf)

compare_vectors_pairwise <- function(list1, list2) {
  # Ensure the lists have the same length
  if (length(list1) != length(list2)) {
    stop("Lists must have the same length.")
  }
  
  # Initialize an empty data frame to store the results
  comparison_results <- data.frame(Pair = character(0), Shared = integer(0),
                                   UniqueInList1 = integer(0), UniqueInList2 = integer(0))
  
  # Iterate through each pair of vectors and calculate the
  # comparison results
  for (i in 1:length(list1)) {
    vector1 <- list1[[i]]
    vector2 <- list2[[i]]
    
    shared_count <- length(intersect(vector1, vector2))
    unique_count_list1 <- length(setdiff(vector1, vector2))
    unique_count_list2 <- length(setdiff(vector2, vector1))
    
    comparison_results <- rbind(comparison_results, data.frame(Pair = paste("Pair",
                                                                            i), Shared = shared_count, UniqueInList1 = unique_count_list1,
                                                               UniqueInList2 = unique_count_list2))
  }
  
  return(comparison_results)
}

list1 <- list(Adult_TF[1, -1], Adult_TF[3, -1], Adult_TF[4, -1], Adult_TF[5,
                                                                          -1], Adult_TF[8, -1], Adult_TF[9, -1])

list2 <- list(fetal_AC, fetal_BC, fetal_Cone, fetal_HC, fetal_RGC, fetal_Rod)

comparison_results <- compare_vectors_pairwise(list1, list2)
print(comparison_results)
comparison_results$Pair <- c("AC", "BC", "Cone", "HC", "RGC", "Rod")
colnames(comparison_results) <- c("Pair","Shared","Adult Specific","Fetal Specific")


Adult_TF <- read_excel("Desktop/Adult_TF.xlsx",sheet = "Major class RSS")
Adult_TF <- extract_before_underscore(colnames(Adult_TF))[-1]
a<-length(intersect(Adult_TF,unique(c(fetal_AC,fetal_BC,fetal_Cone,fetal_RGC,fetal_RGC,fetal_HC))))
b<-length(setdiff(Adult_TF,unique(c(fetal_AC,fetal_BC,fetal_Cone,fetal_RGC,fetal_RGC,fetal_HC))))
c<-length(setdiff(unique(c(fetal_AC,fetal_BC,fetal_Cone,fetal_RGC,fetal_RGC,fetal_HC)),Adult_TF))
comparison_results <- rbind(comparison_results, c("Union",a,b,c))

data = comparison_results %>%
  pivot_longer(cols = -"Pair")
data <- data.frame(data)
data$value <- as.numeric(data$value)
ggplot(data, aes(x = factor(Pair), y = value, fill = name)) + geom_bar(stat = "identity") +
  labs(title = "Compare Overlap and Unique TFs Between Fetal and Adult", x = "Comparison",
       y = "Count") + theme_minimal() + theme(axis.text.x = element_text(angle = 45,
                                                                         hjust = 1))

