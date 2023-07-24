library(ggplot2)
library(ggrepel)  
data <- read.csv("/Users/zhenzuo/Desktop/test.csv")
data$group = ""
data = data[abs(data$logfoldchanges) > 1,]
special_chars_pattern <- "[[:punct:]]"  # This pattern matches any punctuation character

# Check if the string contains any special characters
contains_special_chars <- grepl(special_chars_pattern, data$names)
data = data[!contains_special_chars,]
data[(data$logfoldchanges>0)&(data$RPE_log2FoldChange>0),"group"] = "Cone/RPE macula enriched"
data[(data$logfoldchanges<0)&(data$RPE_log2FoldChange<0),"group"] = "Cone/RPE peripheral enriched"
data[(data$logfoldchanges>0)&(data$RPE_log2FoldChange<0),"group"] = "Cone macula and RPE peripheral enriched"
data[(data$logfoldchanges<0)&(data$RPE_log2FoldChange>0),"group"] = "Cone peripheral and RPE macula enriched"
p <- ggplot(data) + 
  geom_point(aes(RPE_log2FoldChange, logfoldchanges,color = group),size = 0)+
  geom_text_repel(size = 10, aes(RPE_log2FoldChange, logfoldchanges, label = names,color = group),max.overlaps = Inf,show.legend = FALSE)+theme(text = element_text(family = "Arial"))+labs(y= "", x = "")+
   theme_bw()+ theme(legend.text=element_text(size=20))+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE,title=""))+
  scale_fill_discrete(name = "")+
  guides(color = guide_legend(
    override.aes=list(size =18)))+ theme(legend.position = c(0, 0))+ theme(legend.position="top")
ggsave("/Users/zhenzuo/Desktop/test.svg",p,width = 22,height = 10) s
