t_sw <- read.csv("/Users/zhenzuo/Desktop/t_sw.csv")
rownames(t_sw) <- t_sw$Gene
t_sw <- t_sw[,2:5]
t_sw <- t_sw/20
module <- read.csv("/Users/zhenzuo/Desktop/mod_reordered.csv")
rownames(module) <- module$X
t_sw$module <- module[rownames(t_sw),"Module"]
t_sw
t_sw <- t_sw[rowSums(is.na(t_sw)) != 4, ]
write.csv(t_sw,"/Users/zhenzuo/Desktop/t_sw2.csv")
library(dplyr)
# Define the function
t_sw %>%
  # Specify group indicator, column, function
  group_by(module) %>%
  # Calculate the mean of the "Frequency" column for each group
  summarise_at(vars(primed, coupled.on, decoupled, coupled.off), list(Mean = ~mean(., na.rm = TRUE)))


a <- t_sw$coupled.on[t_sw$module==1]

b <- t_sw$coupled.on[t_sw$module==3]

t.test(a,b)

t_sw_2 <- t_sw[t_sw$module==2,]
